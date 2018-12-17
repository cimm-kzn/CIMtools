# -*- coding: utf-8 -*-
#
#  Copyright 2018 Ramil Nugmanov <stsouko@live.ru>
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
from numpy import array, empty
from pandas import DataFrame, Series, concat
from sklearn.utils import column_or_1d
from ..utils import iter2array, nested_iter_to_2d_array


def make_filtered_transformer(transformer, x_is_2d=False):
    class Transformer:
        def fit(self, x, y, **fit_params):
            x, not_nones, _x = self.__check_x(x)
            y = column_or_1d(y, warn=True)
            return transformer.fit(_x, y[not_nones], **fit_params)

        def fit_transform(self, x, y, **fit_params):
            x, not_nones, _x = self.__check_x(x)
            y = column_or_1d(y, warn=True)
            _x = transformer.fit(_x, y[not_nones], **fit_params).transform(_x)
            return self.__prepare_x(x, _x, not_nones)

        def transform(self, x):
            x, not_nones, _x = self.__check_x(x)
            _x = transformer.transform(_x)
            return self.__prepare_x(x, _x, not_nones)

        @staticmethod
        def __check_x(x):
            if x_is_2d:
                if hasattr(transformer, '_dtype'):
                    x = nested_iter_to_2d_array(x, dtype=transformer._dtype, allow_none=True)
                else:
                    x = nested_iter_to_2d_array(x, allow_none=True)

                if isinstance(x, DataFrame):
                    not_nones = x.isna().mean(axis=1) != 1
                    _x = x.loc[not_nones]
                else:
                    not_nones = (x == array(None)).mean(axis=1) != 1
                    _x = x[not_nones]
            else:
                if hasattr(transformer, '_dtype'):
                    x = iter2array(x, dtype=transformer._dtype, allow_none=True)
                else:
                    x = iter2array(x, allow_none=True)

                if isinstance(x, Series):
                    not_nones = ~x.isna()
                    _x = x.loc[not_nones]
                else:
                    not_nones = x != array(None)
                    _x = x[not_nones]
            return x, not_nones, _x

        @staticmethod
        def __prepare_x(x, _x, not_nones):
            if len(x) == len(_x):
                return _x

            if isinstance(_x, Series):
                out = Series(index=range(len(x)), dtype=_x.dtype, name=_x.name)
                out.loc[not_nones] = _x.values
            elif isinstance(_x, DataFrame):
                out = [Series(index=range(len(x)), dtype=t, name=c) for c, t in zip(_x.columns, _x.dtypes)]
                for s, c in zip(out, _x.columns):
                    s.loc[not_nones] = _x[c].values
                out = concat(out, axis=1)
            else:
                if len(_x.shape) > 1:
                    shape = (len(x), _x.shape[1])
                else:
                    shape = len(x)

                out = empty(shape, dtype=_x.dtype)
                out[not_nones] = _x
                if _x.dtype.kind != 'O':
                    out[~not_nones] = None
            return out

        def __getattribute__(self, name):
            if name in ('fit', 'fit_transform', 'transform', '_Transformer__check_x', '_Transformer__prepare_x'):
                return super().__getattribute__(name)
            return getattr(transformer, name)

        def __repr__(self):
            return repr(transformer)

    return Transformer()


__all__ = ['make_filtered_transformer']
