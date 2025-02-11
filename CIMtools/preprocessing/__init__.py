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
from .cgr import *
from .conditions_container import *
from .equation import *
from .fingerprint import *
from .fragmentor import *
from .graph_encoder import *
from .graph_to_matrix import *
from .solvent import *
from .standardize import *
from .standardize import __all__ as _standardize


__all__ = ['Conditions', 'DictToConditions', 'ConditionsToDataFrame', 'SolventVectorizer', 'EquationTransformer',
           'CGR', 'MoleculesToMatrix', 'CGRToMatrix']
__all__.extend(_standardize)

if 'Fragmentor' in locals():
    __all__.append('Fragmentor')
    __all__.append('FragmentorFingerprint')
if 'GNNFingerprint' in locals():
    __all__.append('GNNFingerprint')
