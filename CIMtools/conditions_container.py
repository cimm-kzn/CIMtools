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
from operator import itemgetter
from .metric_constants import C, bar


class Conditions:
    def __init__(self, temperature=25 * C, pressure=1 * bar, solvents=None):
        self.temperature = temperature
        self.pressure = pressure
        self.solvents = solvents

    @property
    def temperature(self):
        return self.__temperature

    @temperature.setter
    def temperature(self, value):
        value = float(value)
        if value <= 0:
            raise ValueError('only positive temperature possible in this universe')
        self.__temperature = value

    @property
    def pressure(self):
        return self.__pressure

    @pressure.setter
    def pressure(self, value):
        value = float(value)
        if value <= 0:
            raise ValueError('only positive pressure possible in this universe')
        self.__pressure = value

    @property
    def solvents(self):
        return self.__solvents

    @solvents.setter
    def solvents(self, value):
        solvents = []
        for k, v in value:
            if k not in known_solvents:
                raise KeyError('unknown solvent')
            v = float(v)
            if v <= 0 or v > 1:
                raise ValueError('impossible solvent amount')
            solvents.append((k, v))
        if sum(x for _, x in solvents) > 1:
            raise ValueError('impossible total amount of solvents')
        self.__solvents = tuple(sorted(solvents, key=itemgetter(1), reverse=True))


known_solvents = dict((  # DON'T TOUCH ORDER. ONLY APPEND.
    ('1,3-xylene', 'CC1=CC(C)=CC=C1'),
    ('anisole', 'COC1=CC=CC=C1'),
    ('1,2-dichloroethane', 'ClCCCl'),
    ('benzene', 'C1=CC=CC=C1'),
    ('formamide', 'NC=O'),
    ('ethoxyethane', 'CCOCC'),
    ('ethyl benzoate', 'CCOC(=O)C1=CC=CC=C1'),
    ('heptan-1-ol', 'CCCCCCCO'),
    ('3-methylbutan-1-ol', 'CC(C)CCO'),
    ('bromobenzene', 'BrC1=CC=CC=C1'),
    ('butan-1-ol', 'CCCCO'),
    ('butan-2-ol', 'CCC(C)O'),
    ('pentan-1-ol', 'CCCCCO'),
    ('propan-2-ol', 'CC(C)O'),
    ('cyclohexane', 'C1CCCCC1'),
    ('propan-2-one', 'CC(C)=O'),
    ('chlorobenzene', 'ClC1=CC=CC=C1'),
    ('oxolane', 'C1CCOC1'),
    ('methanedithione', 'S=C=S'),
    ('pyridine', 'C1=CC=NC=C1'),
    ('acetic acid', 'CC(O)=O'),
    ('butan-2-one', 'CCC(C)=O'),
    ('acetonitrile', 'CC##N'),
    ('2,2,4-trimethylpentane', 'CC(C)CC(C)(C)C'),
    ('2-methylpropan-1-ol', 'CC(C)CO'),
    ('2-methylpropan-2-ol', 'CC(C)(C)O'),
    ('methanol', 'CO'),
    ('propan-1-ol', 'CCCO'),
    ('methanesulfinylmethane', 'CS(C)=O'),
    ('dichloromethane', 'ClCCl'),
    ('ethanol', 'CCO'),
    ('water', 'O'),
    ('2-methylbutan-2-ol', 'CCC(C)(C)O'),
    ('butanoate', 'CCCC([O-])=O'),
    ('1-phenylethan-1-one', 'CC(=O)C1=CC=CC=C1'),
    ('1,3,5‐trimethylbenzene', 'CC1=CC(C)=CC(C)=C1'),
    ('ethyl acetate', 'CCOC(C)=O'),
    ('1,4-dioxane', 'C1COCCO1'),
    ('tetrachloromethane', 'ClC(Cl)(Cl)Cl'),
    ('ethane-1,2-diol', 'OCCO'),
    ('trichloromethane', 'ClC(Cl)Cl'),
    ('propanenitrile', 'CCC##N'),
    ('phenylmethanol', 'OCC1=CC=CC=C1'),
    ('oxolan-2-one', 'O=C1CCCO1'),
    ('hexan-1-ol', 'CCCCCCO'),
    ('hexane', 'CCCCCC'),
    ('octan-1-ol', 'CCCCCCCCO'),
    ('nitromethane', 'C[N+]([O-])=O'),
    ('nitrobenzene', '[O-][N+](=O)C1=CC=CC=C1'),
    ('toluene', 'CC1=CC=CC=C1'),
    ('heavy water', '[2H]O[2H]'),
    ('tetrahydrothiophene 1,1-dioxide', 'O=S1(=O)CCCC1'),
    ('hexamethylphosphoramide', 'CN(C)P(=O)(N(C)C)N(C)C'),
    ('n,n-dimethylformamide', 'CN(C)C=O'),
    ('n,n-dimethylacetamide', 'CN(C)C(C)=O'),
    ('1,2-dimethoxyethane', 'COCCOC'),
    ('1,2,3-propanetriol', 'OCC(O)CO')))


__all__ = ['Conditions']
