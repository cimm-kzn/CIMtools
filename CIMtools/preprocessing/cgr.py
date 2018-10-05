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
from CGRtools import CGRpreparer
from CGRtools.containers import ReactionContainer
from sklearn.base import BaseEstimator
from .common import iter2array, TransformerMixin
from ..exceptions import ConfigurationError


class CGR(BaseEstimator, TransformerMixin):
    def __init__(self, cgr_type='0', extralabels=False):
        self.cgr_type = cgr_type
        self.extralabels = extralabels
        self.__init()

    def __init(self):
        try:
            self.__cgr = CGRpreparer(cgr_type=self.cgr_type, extralabels=self.extralabels)
        except Exception as e:
            raise ConfigurationError(e)

    def __getstate__(self):
        return {k: v for k, v in super().__getstate__().items() if not k.startswith('_CGR__')}

    def __setstate__(self, state):
        super().__setstate__(state)
        self.__init()

    def set_params(self, **params):
        if params:
            super().set_params(**params)
            self.__init()
        return self

    def transform(self, x):
        x = super().transform(x)

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

    _dtype = ReactionContainer
