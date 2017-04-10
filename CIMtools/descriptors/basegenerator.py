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
from abc import ABC, abstractmethod
from functools import reduce
from operator import and_
from pandas import Series, Index, MultiIndex, concat
from .descriptoragregator import PropertyExtractor


class BaseGenerator(ABC, PropertyExtractor):
    def __init__(self, s_option=None):
        PropertyExtractor.__init__(self, s_option)

    @abstractmethod
    def get_config(self):
        pass

    @abstractmethod
    def prepare(self, structures, **_):
        pass

    @property
    @abstractmethod
    def markers(self) -> int:
        pass

    @abstractmethod
    def set_work_path(self, workpath):
        pass

    @abstractmethod
    def delete_work_path(self):
        pass

    def get(self, structures, **kwargs):
        tmp = self.prepare(structures, **kwargs)
        if not tmp:
            return False

        x, y, ad, i, s = tmp

        res = dict(X=concat(x, axis=1, keys=range(len(x))) if len(x) > 1 else x[0],
                   AD=reduce(and_, ad), Y=Series(y, name='Property'), structures=s)  # todo: prepare structures

        if self.markers:
            _i = MultiIndex.from_tuples(i, names=['structure'] + ['c.%d' % x for x in range(self.markers)])
        else:
            _i = Index(i, name='structure')

        res['X'].index = res['AD'].index = res['Y'].index = _i
        return res

    def write_prepared(self, structures, writers):
        prop = []
        doubles = []
        used_str = []
        for s_numb, s in enumerate(structures):
            if isinstance(s, list):
                meta = s[0][0][1].meta
                for d in s:  # d = ((n1, tmp1), (n2, tmp2), ...)
                    tmp_d = [s_numb]
                    tmp_s = []  # list of graphs with marked atoms
                    for w, (x, y) in zip(writers, d):
                        w.write(y)
                        tmp_d.append(x)
                        tmp_s.append(y)
                    prop.append(self.get_property(meta, marks=tmp_d[1:]))
                    doubles.append(tmp_d)
                    used_str.append(tmp_s)
            else:
                writers[0].write(s)
                prop.append(self.get_property(s.meta))
                doubles.append(s_numb)
                used_str.append(s)

        return prop, doubles, used_str
