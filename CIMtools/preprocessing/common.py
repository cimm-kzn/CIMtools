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
from sklearn.base import TransformerMixin as _TransformerMixin


def iter2array(data, dtype=BaseContainer, allow_none=False):
    if isinstance(data, ndarray):
        assert len(data.shape) == 1, 'invalid input array shape'
    elif not isinstance(data, (list, tuple)):  # try to unpack iterable
        data = list(data)

    assert len(data), 'empty input array'
    if allow_none:
        assert all(isinstance(x, dtype) for x in data if x is not None), 'invalid dtype'
    else:
        assert all(isinstance(x, dtype) for x in data), 'invalid dtype'

    if isinstance(data, ndarray):
        return data

    out = empty(len(data), dtype=object)
    for n, x in enumerate(data):
        if x is not None:
            out[n] = x
    return out


def nested_iter_to_2d_array(data, dtype=BaseContainer, allow_none=False):
    if isinstance(data, ndarray):
        assert len(data.shape) == 2, 'invalid input array shape'
        assert data.size, 'empty input array'
        if not allow_none:
            assert all(isinstance(x, dtype) for x in data.flat), 'invalid dtype'
        return data

    if not isinstance(data, (list, tuple)):  # try to unpack iterable
        data = list(data)
    if not all(isinstance(x, (list, tuple)) for x in data):
        data = [x if isinstance(x, (list, tuple)) else list(x) for x in data]

    assert len(data), 'empty input array'

    if allow_none:
        shape = max((len(x) for x in data), default=0)
        assert shape, 'empty input array'
        assert all(isinstance(y, dtype) for x in data for y in x if y is not None), 'invalid dtype'
    else:
        shape = len(data[0])
        assert shape, 'empty input array'
        assert all(len(x) == shape for x in data[1:]), 'input data contains None'
        assert all(isinstance(y, dtype) for x in data for y in x), 'invalid dtype'

    out = empty((len(data), shape), dtype=object)
    for n, x in enumerate(data):
        for m, y in enumerate(x):
            if y is not None:
                out[n, m] = y
    return out


class TransformerMixin(_TransformerMixin):
    def fit(self, x, y=None):
        """Do nothing and return the estimator unchanged

        This method is just there to implement the usual API and hence work in pipelines.
        """
        iter2array(x)
        return self

    def transform(self, x, y=None):
        return iter2array(x)
