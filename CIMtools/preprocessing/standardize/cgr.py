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
from CGRtools.containers import MoleculeContainer, ReactionContainer
from pandas import DataFrame
from ...base import CIMtoolsTransformerMixin


class StandardizeCGR(CIMtoolsTransformerMixin):
    def __init__(self):
        """
        Reactions and Molecules standardization

        For molecules kekule/thiele and groups standardization procedures will be applied.
        """

    def transform(self, x):
        return DataFrame([[self.__prepare(g)] for g in super().transform(x)], columns=['standardized'])

    @staticmethod
    def __prepare(r):
        r = r.copy()
        r.canonicalize()
        return r

    _dtype = (MoleculeContainer, ReactionContainer)


__all__ = ['StandardizeCGR']
