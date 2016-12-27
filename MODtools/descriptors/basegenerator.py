#!/usr/bin/env python3.4
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
import operator
import pandas as pd
from functools import reduce
from .descriptoragregator import Propertyextractor


class BaseGenerator(Propertyextractor):
    def __init__(self, workpath='.', s_option=None):
        Propertyextractor.__init__(self, s_option)

        self.workpath = workpath

    def setworkpath(self, workpath):
        self.workpath = workpath

    def get(self, structures, **kwargs):
        tmp = self.prepare(structures, **kwargs)
        if not tmp:
            return False

        X, Y, AD, I, S = tmp

        res = dict(X=pd.concat(X, axis=1, keys=range(len(X))) if len(X) > 1 else X[0],
                   AD=reduce(operator.and_, AD), Y=pd.Series(Y, name='Property'),
                   structures=S)  # todo: prepare structures

        if self.markers:
            i = pd.MultiIndex.from_tuples(I, names=['structure'] + ['c.%d' % x for x in range(self.markers)])
        else:
            i = pd.Index(I, name='structure')

        res['X'].index = res['AD'].index = res['Y'].index = i
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
