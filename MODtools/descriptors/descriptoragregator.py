# -*- coding: utf-8 -*-
#
# Copyright 2016 Ramil Nugmanov <stsouko@live.ru>
# This file is part of MODtools.
#
# MODtools is free software; you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Affero General Public License for more details.
#
#  You should have received a copy of the GNU Affero General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#
import sys
import operator
import pandas as pd
from collections import defaultdict
from functools import reduce
from CGRtools.files.SDFrw import SDFread
from CGRtools.files.RDFrw import RDFread


class Descriptorchain(object):
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

    def setworkpath(self, workpath):
        for gen, _ in self.__generators:
            if hasattr(gen, 'setworkpath'):
                gen.setworkpath(workpath)

    def get(self, structures, **kwargs):
        """
        :param structures: opened structure file or stringio
        :param kwargs: generators specific arguments
        :return: dict(X=DataFrame, AD=Series, Y=Series, BOX=Series, structures=DataFrame)
        """
        res = defaultdict(list)

        def merge_wrap(x, y):
            return pd.merge(x, y, how='outer', left_index=True, right_index=True)

        for gen, ad in self.__generators:
            for k, v in gen.get(structures, **kwargs).items():
                res[k].append(v)
            res['BOX'].append(pd.Series(ad, index=res['X'][-1].columns))
            structures.seek(0)

        res['X'] = reduce(merge_wrap, res['X'])
        res['AD'] = reduce(operator.and_, sorted(res['AD'], key=lambda x: len(x.index), reverse=True))
        # на данный момент не придумано как поступать с мультицентровостью. пока свойство просто дублируется.
        res['Y'] = sorted(res['Y'], key=lambda x: len(x.index), reverse=True)[0]

        res['BOX'] = pd.concat(res['BOX'])

        if 'structures' in res:
            res['structures'] = reduce(merge_wrap, res['structures'])

        return dict(res)


class Propertyextractor(object):
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


class Descriptorsdict(Propertyextractor):
    def __init__(self, data=None, s_option=None, is_reaction=False, ):
        Propertyextractor.__init__(self, s_option)
        self.__is_reaction = is_reaction
        self.__extention = data
        self.__extheader = self.__prepareextheader(data)

    @staticmethod
    def setworkpath(_):
        return None

    @staticmethod
    def __prepareextheader(data):
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

    def __parsefile(self, structures):
        """
        parse SDF or RDF on known keys-headers.
        :param structures: opened file
        :return: DataFrame of descriptors. indexes is the numbers of structures in file, columns - names of descriptors
        """
        extblock = []
        props = []
        reader = RDFread(structures) if self.__is_reaction else SDFread(structures)
        for i in reader.read():
            meta = i['meta'] if self.__is_reaction else i.graph['meta']
            tmp = []
            for key, value in meta.items():
                if key in self.__extention:
                    data = self.__extention[key]['value'].loc[self.__extention[key]['key'] == value] if \
                        self.__extention[key] else pd.DataFrame([{key: float(value)}])
                    if not data.empty:
                        data.index = [0]
                        tmp.append(data)
            extblock.append(pd.concat(tmp, axis=1) if tmp else pd.DataFrame([{}]))

            props.append(self.get_property(meta))

        res = pd.DataFrame(pd.concat(extblock), columns=self.__extheader)
        res.index = pd.Index(range(len(res.index)), name='structure')
        prop = pd.Series(props, name='Property', index=res.index)
        return res, prop

    def __parseadditions0(self, **kwargs):
        extblock = []
        for i, j in kwargs.items():
            if i in self.__extention:
                for n, k in enumerate(j) if isinstance(j, list) else j.items():
                    data = self.__extention[i]['value'].loc[self.__extention[i]['key'] == k] if \
                        self.__extention[i] else pd.DataFrame([{i: k}])
                    if not data.empty:
                        data.index = [0]
                        if len(extblock) > n:
                            extblock[n].append(data)
                        else:
                            extblock.extend([[] for _ in range(n - len(extblock))] + [[data]])
        res = pd.DataFrame(pd.concat([pd.concat(x, axis=1) if x else pd.DataFrame([{}]) for x in extblock]),
                           columns=self.__extheader)
        res.index = pd.Index(range(len(res.index)), name='structure')
        return res

    def __parseadditions1(self, **kwargs):
        tmp = []
        for i, j in kwargs.items():
            if i in self.__extention:
                data = self.__extention[i]['value'].loc[self.__extention[i]['key'] == j] if \
                       self.__extention[i] else pd.DataFrame([{i: j}])
                if not data.empty:
                    data.index = [0]
                    tmp.append(data)
        return pd.DataFrame(pd.concat(tmp, axis=1) if tmp else pd.DataFrame([{}]), columns=self.__extheader,
                            index=pd.Index([0], name='structure'))

    def get(self, structures=None, **kwargs):
        if kwargs.get('parsesdf'):
            extblock, prop = self.__parsefile(structures)
            structures.seek(0)  # ad-hoc for rereading

        elif all(isinstance(x, list) or isinstance(x, dict) for y, x in kwargs.items() if y in self.__extention):
            extblock = self.__parseadditions0(**kwargs)
            prop = pd.Series(index=extblock.index)

        elif not any(isinstance(x, list) or isinstance(x, dict) for y, x in kwargs.items() if y in self.__extention):
            extblock = self.__parseadditions1(**kwargs)
            prop = pd.Series(index=extblock.index)

        else:
            print('WHAT DO YOU WANT? use correct extentions params', file=sys.stderr)
            return False

        return dict(X=extblock, AD=-extblock.isnull().any(axis=1), Y=prop)
