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
from collections import defaultdict
from .score import MultiLabelScore


def iap(y_test, y_prob):
    res = defaultdict(list)
    for sid, s_prob in y_prob.groupby(level='structure'):
        s_test = y_test.xs(sid, level='structure', drop_level=False)
        for col in s_prob.columns:
            in_class = s_test.loc[s_test == col].index
            out_class = s_test.loc[s_test != col].index

            in_class_dominance = sum(s_prob.loc[i][col] > s_prob.loc[o][col] for i, o in product(in_class, out_class))
            possible_combo = len(in_class) * len(out_class)
            if possible_combo:
                res[col].append(in_class_dominance / possible_combo)
    res = MultiLabelScore({x: sum(y) / len(y) for x, y in res.items()})
    return res
