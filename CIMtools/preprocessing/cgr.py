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
from CGRtools import CGRpreparer
from CGRtools.containers import ReactionContainer
from sklearn.base import BaseEstimator
from .common import iter2array, TransformerMixin


class CGR(BaseEstimator, TransformerMixin):
    def __init__(self, cgr_type='0', extralabels=False):
        self.cgr_type = cgr_type
        self.extralabels = extralabels

    def __getstate__(self):
        return {k: v for k, v in super().__getstate__().items() if not k.startswith('_CGR__')}

    def set_params(self, **params):
        if params:
            super().set_params(**params)
            self.__cgr = None
        return self

    def transform(self, x):
        x = super().transform(x)

        if self.__cgr is None:
            self.__cgr = CGRpreparer(cgr_type=self.cgr_type, extralabels=self.extralabels)

        res = []
        for n, s in enumerate(x):
            try:
                c = self.__cgr.condense(s)
            except Exception as e:
                print('Skip reaction:', n, e)
                res.append(None)
            else:
                res.append(c)

        return iter2array(res, allow_none=True)

    __cgr = None
    _dtype = ReactionContainer
