# -*- coding: utf-8 -*-
#
#  Copyright 2020 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from CGRtools.containers import ReactionContainer, MoleculeContainer
from pathlib import Path
from standardize import StandardizeChemAxon
from . import __path__
from ..exceptions import ConfigurationError
from ..utils import iter2array


class MappingChemAxon(StandardizeChemAxon):
    def __init__(self):
        pass
 
    def transform(self, x):
        x = iter2array(x, dtype=(MoleculeContainer, ReactionContainer))
        
        rules = open(Path(__path__[0]) /'map_reaction.xml').read()
        chem_axon_rule = StandardizeChemAxon(rules)
        
        try:
            mapped = chem_axon_rule.fit_transform(x)
        except Exception as e:
            raise ConfigurationError(e)
        else:
            if len(mapped) != len(x):
                raise ValueError('invalid data')
        return mapped

__all__ = ['MappingChemAxon']
