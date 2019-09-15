# -*- coding: utf-8 -*-
#
#  Copyright 2018, 2019 Ramil Nugmanov <stsouko@live.ru>
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
from CGRtools.reactor import CGRReactor
from CGRtools.containers import MoleculeContainer, CGRContainer, ReactionContainer
from sklearn.base import BaseEstimator
from ...base import CIMtoolsTransformerMixin
from ...exceptions import ConfigurationError
from ...utils import iter2array


class StandardizeCGR(BaseEstimator, CIMtoolsTransformerMixin):
    def __init__(self, templates=(), delete_atoms=False):
        """
        Molecule and CGR standardization

        For molecules kekule/thiele and groups standardization procedures will be applied.

        :param templates: CGRTemplates. list of rules for graph modifications.
        :param delete_atoms: if True atoms exists in templates reactants but not exists in products will be removed
        """
        self.templates = templates
        self.delete_atoms = delete_atoms
        self.__init()

    def __init(self):
        try:
            self.__fixes = [CGRReactor(x, delete_atoms=self.delete_atoms) for x in self.templates]
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
        return iter2array(self.__prepare(g) for g in super().transform(x))

    def __prepare(self, g):
        if isinstance(g, MoleculeContainer):
            g = g.copy()
            g.standardize()
            g.kekule()
            g.thiele()

        for fix in self.__fixes:
            while True:
                try:
                    p = next(fix(g, False))
                except StopIteration:
                    break
                else:
                    p.meta.update(g.meta)
                    g = p
        return g

    _dtype = (MoleculeContainer, CGRContainer)


__all__ = ['StandardizeCGR']
