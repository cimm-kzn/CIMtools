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
from CGRtools import CGRreactor
from CGRtools.containers import MoleculeContainer
from sklearn.base import BaseEstimator
from ...base import CIMtoolsTransformerMixin
from ...exceptions import ConfigurationError
from ...utils import iter2array


class StandardizeCGR(BaseEstimator, CIMtoolsTransformerMixin):
    def __init__(self, templates, balance_groups=False):
        """
        CGR standardization and reaction balancing

        :param templates: CGRTemplates. list of rules for graph modifications.
        :param balance_groups: if True: for unbalanced reactions contains multiple attached functional groups in
            products and one of them described in reagents - will be restored information about all equal groups.
            for example:

                R + B1-X-> B'1-R'-B'2 + X'

            where B' is transformed B, R and X same.
            we know what B'1 and B'2 is equal and B'1 is transformed B1 =>
            this groups most likely appeared from a single reagent. we can add copy of B-X to reagents.
            results will be:

                R + B1-X1 + B2-X2 -> B'1-R'-B'2 + X'1 + X'2

        """
        self.templates = templates
        self.balance_groups = balance_groups
        self.__init()

    def __init(self):
        try:
            assert self.templates or self.balance_groups, 'invalid params. need balance_groups or/and templates'
            self.__fixes = [CGRreactor(x) for x in self.templates]
        except Exception as e:
            raise ConfigurationError from e

    def __getstate__(self):
        return {k: v for k, v in super().__getstate__().items() if not k.startswith('_StandardizeCGR__')}

    def __setstate__(self, state):
        super().__setstate__(state)
        self.__init()

    def set_params(self, **params):
        if params:
            super().set_params(**params)
            self.__init()
        return self

    def transform(self, x):
        return iter2array((self.__prepare(g) for g in super().transform(x)), allow_none=True)

    def __prepare(self, g):
        if self.balance_groups:  # DO NOT WORKING
            g = clone_subgraphs(g)

        for fix in self.__fixes:
            while True:
                p = fix(g)
                if p:
                    p.meta.update(g.meta)
                    g = p
                else:
                    break
        return g

    _dtype = MoleculeContainer


__all__ = ['StandardizeCGR']
