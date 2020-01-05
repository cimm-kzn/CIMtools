# -*- coding: utf-8 -*-
#
#  Copyright 2018, 2019 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from CGRtools import CGRPreparer
from CGRtools.containers import ReactionContainer
from sklearn.base import BaseEstimator
from ..base import CIMtoolsTransformerMixin
from ..exceptions import ConfigurationError
from ..utils import iter2array


class CGR(BaseEstimator, CIMtoolsTransformerMixin):
    def __init__(self, cgr_type='0'):
        self.cgr_type = cgr_type
        self.__init()

    def __init(self):
        try:
            self.__cgr = CGRPreparer(self.cgr_type)
        except Exception as e:
            raise ConfigurationError from e

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
        cgr = self.__cgr
        return iter2array(cgr.compose(s) for s in x)

    _dtype = ReactionContainer


__all__ = ['CGR']
