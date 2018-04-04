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
from numpy import array, empty
from sklearn.utils import column_or_1d
from .preprocessing.common import iter2array, nested_iter_to_2d_array


def make_filtered_transformer(transformer, x_is_2d=False):
    class Transformer:
        def fit(self, x, y=None, **fit_params):
            x, not_nones = self.__check_x(x)
            y = column_or_1d(y, warn=True)
            return transformer.fit(x[not_nones], y[not_nones], **fit_params)

        def fit_transform(self, x, y=None, **fit_params):
            x, not_nones = self.__check_x(x)
            y = column_or_1d(y, warn=True)
            _x = transformer.fit(x[not_nones], y[not_nones], **fit_params).transform(x[not_nones])
            return self.__prepare_x(x, _x, not_nones)

        def transform(self, x):
            x, not_nones = self.__check_x(x)
            _x = transformer.transform(x[not_nones])
            return self.__prepare_x(x, _x, not_nones)

        @staticmethod
        def __check_x(x):
            if x_is_2d:
                if hasattr(transformer, '_dtype'):
                    x = nested_iter_to_2d_array(x, dtype=transformer._dtype, allow_none=True)
                else:
                    x = nested_iter_to_2d_array(x, allow_none=True)

                not_nones = (x == array(None)).mean(axis=1) != 1
            else:
                if hasattr(transformer, '_dtype'):
                    x = iter2array(x, dtype=transformer._dtype, allow_none=True)
                else:
                    x = iter2array(x, allow_none=True)

                not_nones = x != array(None)
            return x, not_nones

        @staticmethod
        def __prepare_x(x, _x, not_nones):
            if len(_x.shape) > 1:
                shape = (x.shape[0], _x.shape[1])
            else:
                shape = x.shape[0]

            out = empty(shape, dtype=_x.dtype)
            out[not_nones] = _x
            if _x.dtype.kind != 'O':
                out[~not_nones] = None
            return out

        def __getattribute__(self, name):
            if name in ('fit', 'fit_transform', 'transform', '__check_x', '__prepare_x'):
                return super().__getattribute__(name)
            return getattr(transformer, name)

        def __repr__(self):
            return repr(transformer)

    return Transformer()
