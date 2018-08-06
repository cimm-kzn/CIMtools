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
from CGRtools.containers import ReactionContainer, MoleculeContainer
from itertools import tee, chain
from numbers import Number
from numpy import ndarray
from pandas import DataFrame, Series
from sklearn.base import TransformerMixin as _TransformerMixin


def iter2array(data, dtype=(MoleculeContainer, ReactionContainer), allow_none=False):
    if isinstance(data, ndarray):
        assert len(data.shape) == 1, 'invalid input array shape'
    elif not isinstance(data, (Series, list, tuple)):  # try to unpack iterable
        data = list(data)

    assert len(data), 'empty input array'
    if allow_none:
        assert all(isinstance(x, dtype) for x in data if x is not None), 'invalid dtype'
    else:
        assert all(isinstance(x, dtype) for x in data), 'invalid dtype'

    if isinstance(data, (ndarray, Series)):
        return data

    if isinstance(dtype, tuple):
        dtype = all(issubclass(x, Number) for x in dtype) and dtype[0] or object
    if not issubclass(dtype, Number):
        dtype = object

    return Series(data, dtype=dtype)


def nested_iter_to_2d_array(data, dtype=(MoleculeContainer, ReactionContainer), allow_none=False):
    if isinstance(data, (ndarray, DataFrame)):
        assert data.size, 'empty input array'
        if isinstance(data, DataFrame):
            flat = data.values.flat
        else:
            assert len(data.shape) == 2, 'invalid input array shape'
            flat = data.flat

        if allow_none:
            assert all(isinstance(x, dtype) for x in flat if x is not None), 'invalid dtype'
        else:
            assert all(isinstance(x, dtype) for x in flat), 'invalid dtype'
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

    if isinstance(dtype, tuple):
        dtype = all(issubclass(x, Number) for x in dtype) and dtype[0] or object
    if not issubclass(dtype, Number):
        dtype = object

    return DataFrame(data, dtype=dtype)


class TransformerMixin(_TransformerMixin):
    def fit(self, x, y=None):
        """Do nothing and return the estimator unchanged

        This method is just there to implement the usual API and hence work in pipelines.
        """
        if self._dtype is not None:
            iter2array(x, dtype=self._dtype)
        else:
            iter2array(x)
        return self

    def transform(self, x):
        if self._dtype is not None:
            return iter2array(x, dtype=self._dtype)
        return iter2array(x)

    _dtype = None


def reaction_support(_class):
    class ReactionSupport(_class):
        def transform(self, x, y=None):
            assert all(isinstance(s, ReactionContainer) for s in x), 'invalid dtype, olny ReactionContainers acceptable'

            shifts = {}
            mols = []
            for i in ('reagents', 'products'):
                sh = shifts[i] = [len(mols)]
                for s in x:
                    si = s[i]
                    sh.append(len(si) + sh[-1])
                    mols.extend(si)

            transformed = super().transform(mols, y=y)
            assert len(transformed) == len(mols), 'unexpected transformed molecules amount'

            out = []
            for s, (r, p) in zip(x, (transformed[y: z] for y, z in self.__pairwise(shifts['reagents'])),
                                 (transformed[y: z] for y, z in self.__pairwise(shifts['products']))):
                if any(i is None for i in chain(r, p)):
                    out.append(None)
                else:
                    out.append(ReactionContainer(r, p, meta=s.meta))
            return iter2array(out, allow_none=True)

        @staticmethod
        def __pairwise(iterable):
            """s -> (s0,s1), (s1,s2), (s2, s3), ..."""
            a, b = tee(iterable)
            next(b, None)
            return zip(a, b)

    return ReactionSupport
