# -*- coding: utf-8 -*-
#
#  Copyright 2016-2018 Ramil Nugmanov <stsouko@live.ru>
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
from pandas import concat


def save_svm(outputfile, x, y, header=True):
    with open(outputfile + '.svm', 'w', encoding='utf-8') as f:
        if header:
            f.write(' '.join(['Property'] + ['%s:%s' % i for i in enumerate(x.columns, start=1)]) + '\n')

        for i, j in zip(x.values, y):
            f.write(' '.join(['%s ' % j] + ['%s:%s' % x for x in enumerate(i, start=1) if x[1] != 0]) + '\n')


def save_csv(outputfile, x, y, header=True):
    concat([y, x], axis=1).to_csv(outputfile + '.csv', index=False, header=header)
