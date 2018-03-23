# -*- coding: utf-8 -*-
#
#  Copyright 2018 Ramil Nugmanov <stsouko@live.ru>
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
from CGRtools.containers.common import BaseContainer
from numpy import empty, ndarray


def iter2array(data, dtype=BaseContainer):
    if isinstance(data, ndarray):
        assert len(data.shape) == 1, 'invalid input array shape'
    elif not isinstance(data, (list, set, tuple)):
        data = list(data)

    assert data, 'empty input array'
    assert all(isinstance(x, dtype) for x in data), ''

    if isinstance(data, ndarray):
        return data

    out = empty(len(data), dtype=object)
    for n, x in enumerate(data):
        out[n] = x

    return out
