# -*- coding: utf-8 -*-
#
#  Copyright 2016, 2017 Ramil Nugmanov <stsouko@live.ru>
#  This file is part of CIMtools.
#
#  CIMtools is free software; you can redistribute it and/or modify
#  it under the terms of the GNU Affero General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Affero General Public License for more details.
#
#  You should have received a copy of the GNU Affero General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#
from sys import stderr
from operator import and_
from pandas import DataFrame, Series, concat, merge, Index
from collections import defaultdict
from functools import reduce


class DescriptorsChain(object):
    def __init__(self, *args):
        """
        chaining multiple descriptor generators.
        concatenate X vectors and merge AD
        :param args: set of generators or set of list[generator, consider their AD {True|False}]
        """
        if isinstance(args[0], tuple):
            self.__generators = args
        else:
            self.__generators = [(x, True) for x in args]

    def set_work_path(self, workpath):
        for gen, _ in self.__generators:
            if hasattr(gen, 'set_work_path'):
                gen.set_work_path(workpath)

    def flush(self):
        for gen, _ in self.__generators:
            if hasattr(gen, 'flush'):
                gen.flush()

    def get(self, structures, **kwargs):
        """
        :param structures: list of CGRtools data
        :param kwargs: generators specific arguments
        :return: dict(X=DataFrame, AD=Series, Y=Series, BOX=Series, structures=DataFrame)
        """
        res = defaultdict(list)

        def merge_wrap(x, y):
            return merge(x, y, how='outer', left_index=True, right_index=True)

        for gen, ad in self.__generators:
            for k, v in gen.get(structures, **kwargs).items():
                res[k].append(v)
            res['BOX'].append(Series(ad, index=res['X'][-1].columns))

        res['X'] = reduce(merge_wrap, res['X'])
        res['AD'] = reduce(and_, sorted(res['AD'], key=lambda x: len(x.index), reverse=True))
        # на данный момент не придумано как поступать с мультицентровостью. пока свойство просто дублируется.
        res['Y'] = sorted(res['Y'], key=lambda x: len(x.index), reverse=True)[0]

        res['BOX'] = concat(res['BOX'])

        if 'structures' in res:
            res['structures'] = reduce(merge_wrap, res['structures'])

        return dict(res)


class PropertyExtractor(object):
    def __init__(self, name):
        self.__name = name

    def get_property(self, meta, marks=None):
        """
        for marked atom property can named property_name.1-2-3 - where 1-2-3 sorted marked atoms.
        for correct work NEED in rdf mapping started from 1 without breaks.
        or used common property with key property_name
        :param meta: dict of data
        :param marks: list of marked atoms
        :return: float property or None
        """
        tmp = marks and meta.get('%s.%s' % (self.__name,
                                            '-'.join(str(x) for x in sorted(marks)))) or meta.get(self.__name)

        return float(tmp) if tmp else None


class DescriptorsDict(PropertyExtractor):
    def __init__(self, data=None, s_option=None):
        PropertyExtractor.__init__(self, s_option)
        self.__extension = data
        self.__ext_header = self.__prepare_ext_header(data)

    @staticmethod
    def __prepare_ext_header(data):
        """
        :param data: dict
        :return: list of strings. descriptors header
        """
        tmp = []
        for i, j in data.items():
            if j:
                tmp.extend(j['value'].columns)
            else:
                tmp.append(i)
        return tmp

    def __parse_file(self, structures):
        """
        parse SDF or RDF on known keys-headers.
        :param structures: opened file
        :return: DataFrame of descriptors. indexes is the numbers of structures in file, columns - names of descriptors
        """
        extblock = []
        props = []
        for i in structures:
            tmp = []
            for key, value in i.meta.items():
                if key in self.__extension:
                    data = self.__extension[key]['value'].loc[self.__extension[key]['key'] == value] if \
                        self.__extension[key] else DataFrame([{key: float(value)}])
                    if not data.empty:
                        data.index = [0]
                        tmp.append(data)
            extblock.append(concat(tmp, axis=1) if tmp else DataFrame([{}]))

            props.append(self.get_property(i.meta))

        res = DataFrame(concat(extblock), columns=self.__ext_header)
        res.index = Index(range(len(res.index)), name='structure')
        prop = Series(props, name='Property', index=res.index)
        return res, prop

    def __parse_additions_multi(self, **kwargs):
        extblock = []
        for i, j in kwargs.items():
            if i in self.__extension:
                for n, k in enumerate(j) if isinstance(j, list) else j.items():
                    data = self.__extension[i]['value'].loc[self.__extension[i]['key'] == k] if \
                        self.__extension[i] else DataFrame([{i: k}])
                    if not data.empty:
                        data.index = [0]
                        if len(extblock) > n:
                            extblock[n].append(data)
                        else:
                            extblock.extend([[] for _ in range(n - len(extblock))] + [[data]])
        res = DataFrame(concat([concat(x, axis=1) if x else DataFrame([{}]) for x in extblock]),
                        columns=self.__ext_header)
        res.index = Index(range(len(res.index)), name='structure')
        return res

    def __parse_additions_single(self, **kwargs):
        tmp = []
        for i, j in kwargs.items():
            if i in self.__extension:
                data = self.__extension[i]['value'].loc[self.__extension[i]['key'] == j] if \
                       self.__extension[i] else DataFrame([{i: j}])
                if not data.empty:
                    data.index = [0]
                    tmp.append(data)
        return DataFrame(concat(tmp, axis=1) if tmp else DataFrame([{}]), columns=self.__ext_header,
                         index=Index([0], name='structure'))

    def get(self, structures=None, **kwargs):
        if kwargs.get('parsesdf'):
            extblock, prop = self.__parse_file(structures)

        elif all(isinstance(x, list) or isinstance(x, dict) for y, x in kwargs.items() if y in self.__extension):
            extblock = self.__parse_additions_multi(**kwargs)
            prop = Series(index=extblock.index)

        elif not any(isinstance(x, list) or isinstance(x, dict) for y, x in kwargs.items() if y in self.__extension):
            extblock = self.__parse_additions_single(**kwargs)
            prop = Series(index=extblock.index)

        else:
            print('WHAT DO YOU WANT? use correct extentions params', file=stderr)
            return False

        return dict(X=extblock, AD=-extblock.isnull().any(axis=1), Y=prop)
