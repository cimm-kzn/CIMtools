# -*- coding: utf-8 -*-
#
#  Copyright 2018-2020 Ramil Nugmanov <nougmanoff@protonmail.com>
#  This file is part of CIMtools.
#
#  CIMtools is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, see <https://www.gnu.org/licenses/>.
#
from CGRtools.containers import ReactionContainer, MoleculeContainer, CGRContainer
from numbers import Number
from numpy import ndarray, ravel
from pandas import DataFrame, Series


def iter2array(data, dtype=(MoleculeContainer, ReactionContainer, CGRContainer)):
    if isinstance(data, ndarray):
        if len(data.shape) != 1:
            if len(data.shape) == 2 and data.shape[1] == 1:
                data = ravel(data)
            else:
                raise ValueError('invalid data shape')
    elif isinstance(data, DataFrame):
        if data.shape[1] != 1:
            raise ValueError('invalid data shape')
        data = ravel(data)
    elif not isinstance(data, (Series, list, tuple)):  # try to unpack iterable
        data = list(data)

    if not len(data):
        raise ValueError('empty input array')
    if not all(isinstance(x, dtype) for x in data):
        raise TypeError('invalid dtype')

    if isinstance(data, (ndarray, Series)):
        return data

    if isinstance(dtype, tuple):
        dtype = all(issubclass(x, Number) for x in dtype) and dtype[0] or object
    if not issubclass(dtype, Number):
        dtype = object

    return Series(data, dtype=dtype)


def nested_iter_to_2d_array(data, dtype=(MoleculeContainer, ReactionContainer, CGRContainer)):
    if isinstance(data, (ndarray, DataFrame)):
        if not data.size:
            raise ValueError('empty input array')
        if isinstance(data, DataFrame):
            flat = data.values.flat
        else:
            if len(data.shape) != 2:
                raise ValueError('invalid input array shape')
            flat = data.flat

        if not all(isinstance(x, dtype) for x in flat):
            raise TypeError('invalid dtype')
        return data

    if not isinstance(data, (list, tuple)):  # try to unpack iterable
        data = list(data)
    if not data:
        raise ValueError('empty input array')
    if not all(isinstance(x, (list, tuple)) for x in data):
        data = [x if isinstance(x, (list, tuple)) else list(x) for x in data]

    shape = len(data[0])
    if not shape:
        raise ValueError('empty input array')
    elif not all(len(x) == shape for x in data[1:]):
        raise ValueError('input data contains None')
    elif not all(isinstance(y, dtype) for x in data for y in x):
        raise TypeError('invalid dtype')

    if isinstance(dtype, tuple):
        dtype = all(issubclass(x, Number) for x in dtype) and dtype[0] or object
    if not issubclass(dtype, Number):
        dtype = object

    return DataFrame(data, dtype=dtype)


__all__ = ['iter2array', 'nested_iter_to_2d_array']
