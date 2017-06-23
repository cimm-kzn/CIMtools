# -*- coding: utf-8 -*-
#
#  Copyright 2017 Ramil Nugmanov <stsouko@live.ru>
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
from pandas import DataFrame, Series, concat, Index
from sys import stderr
from .basegenerator import DataContainer
from .propertyextractor import PropertyExtractor


class DescriptorsDict(PropertyExtractor):
    def __init__(self, data=None, s_option=None):
        PropertyExtractor.__init__(self, s_option)
        self.__extension = data
        self.__ext_header = self.__prepare_ext_header(data)
        self.__config = dict(data=data, s_option=s_option)

    def pickle(self):
        return self.__config.copy()

    @classmethod
    def unpickle(cls, config):
        if {'data', 's_option'}.difference(config):
            raise Exception('Invalid config')
        return cls(data=config['data'], s_option=config['s_option'])

    @staticmethod
    def __prepare_ext_header(data):
        """
        :param data: dict
        :return: list of strings. descriptors header
        """
        tmp = []
        for i in sorted(data):
            j = data[i]
            if j:
                tmp.extend(sorted(list(j.values())[0]))
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
                    data = DataFrame(self.__extension[key][value] if self.__extension[key] else {key: float(value)},
                                     index=[0])
                    if not data.empty:
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
                    data = DataFrame(self.__extension[i][k] if self.__extension[i] else {i: k}, index=[0])
                    if not data.empty:
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
                data = DataFrame(self.__extension[i][j] if self.__extension[i] else {i: j}, index=[0])
                if not data.empty:
                    tmp.append(data)
        return DataFrame(concat(tmp, axis=1) if tmp else DataFrame([{}]), columns=self.__ext_header,
                         index=Index([0], name='structure'))

    def get(self, structures=None, **kwargs):
        if kwargs.get('in_structures'):
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

        return DataContainer(X=extblock, Y=prop, AD=-extblock.isnull().any(axis=1))
