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
from itertools import chain
from operator import itemgetter
from pandas import DataFrame
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.utils import column_or_1d
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


class DictToConditions:
    def __init__(self, temperature=None, pressure=None, solvents=None, amounts=None,
                 default_temperature=25 * C, default_pressure=1 * bar,
                 default_first_solvent='water', default_first_amount=1):
        """Dictionary to Conditions mapper

        :param temperature: name of temperature key
        :param pressure: name of pressure key
        :param solvents: names of solvents keys
        :param amounts: names of solvents amounts keys
        """
        if solvents:
            if 1 < len(solvents) != len(amounts):
                raise ValueError('solvents and amounts lists should be equal')
        elif amounts and len(amounts) > 1:
            raise ValueError('multiple amount without solvents impossible')

        self.__temperature = temperature
        self.__pressure = pressure
        self.__solvents = solvents
        self.__amounts = amounts
        self.__default_temperature = default_temperature
        self.__default_pressure = default_pressure
        self.__default_first_solvent = default_first_solvent
        self.__default_first_amount = default_first_amount

    def transform(self, x):
        if self.__solvents:
            if len(self.__solvents) > 1:
                solvents = self.__solvents
                solvents = [[d[k] for k in solvents] for d in x]
                amounts = self.__amounts
                amounts = [[d[k] for k in amounts] for d in x]
            else:
                solvents = self.__solvents[0]
                solvents = [[d.get(solvents, self.__default_first_solvent)] for d in x]
                if self.__amounts:
                    amounts = self.__amounts[0]
                    amounts = [[d.get(amounts, self.__default_first_amount)] for d in x]
                else:
                    amounts = [[self.__default_first_amount]] * len(x)
        elif self.__amounts:
            amounts = self.__amounts[0]
            amounts = [[d.get(amounts, self.__default_first_amount)] for d in x]
            solvents = [[self.__default_first_solvent]] * len(x)
        else:
            solvents = amounts = [[]] * len(x)

        if self.__temperature:
            temperatures = self.__temperature
            temperatures = [d.get(temperatures, self.__default_temperature) for d in x]
        else:
            temperatures = [self.__default_temperature] * len(x)

        if self.__pressure:
            pressures = self.__pressure
            pressures = [d.get(pressures, self.__default_pressure) for d in x]
        else:
            pressures = [self.__default_pressure] * len(x)

        return [Conditions(t, p, zip(*s)) for t, p, s in zip(temperatures, pressures, zip(solvents, amounts))]


class ConditionsToDataFrame(BaseEstimator, TransformerMixin):
    def __init__(self, max_solvents=1):
        self.max_solvents = max_solvents

    def get_feature_names(self):
        """Get feature names.

        Returns
        -------
        feature_names : list of strings
            Names of the features produced by transform.
        """
        return ['temperature', 'pressure'] + [f'solvent.{x}' for x in range(1, self.max_solvents + 1)] + \
               [f'solvent_amount.{x}' for x in range(1, self.max_solvents + 1)]

    def transform(self, x):
        x = column_or_1d(x, warn=True)
        res = []
        for x in x:
            solvents, amounts = zip(*x.solvents[:self.max_solvents])
            res.append(chain((x.temperature, x.pressure), solvents, amounts))
        return DataFrame(res, columns=self.get_feature_names())

    def fit(self, x, y=None):
        return self


known_solvents = dict((  # DON'T TOUCH ORDER. ONLY APPEND.
    ('1-phenylethan-1-one', 'CC(=O)c1ccccc1'),
    ('1,2-dichloroethane', 'ClCCCl'),
    ('1,2-dimethoxyethane', 'COCCOC'),
    ('1,2,3-propanetriol', 'OCC(O)CO'),
    ('1,3-xylene', 'CC1=CC(C)=CC=C1'),
    ('1,3,5‐trimethylbenzene', 'CC1=CC(C)=CC(C)=C1'),
    ('1,4-dioxane', 'C1COCCO1'),
    ('2-methylbutan-2-ol', 'CCC(C)(C)O'),
    ('2-methylpropan-1-ol', 'CC(C)CO'),
    ('2-methylpropan-2-ol', 'CC(C)(C)O'),
    ('2,2,4-trimethylpentane', 'CC(C)CC(C)(C)C'),
    ('3-methylbutan-1-ol', 'CC(C)CCO'),
    ('acetic acid', 'CC(O)=O'),
    ('acetonitrile', 'CC##N'),
    ('benzene', 'C1=CC=CC=C1'),
    ('benzonitrile', 'N#Cc1ccccc1'),
    ('bromobenzene', 'BrC1=CC=CC=C1'),
    ('butan-1-ol', 'CCCCO'),
    ('butan-2-ol', 'CCC(C)O'),
    ('butan-2-one', 'CCC(C)=O'),
    ('chlorobenzene', 'ClC1=CC=CC=C1'),
    ('cyclohexane', 'C1CCCCC1'),
    ('dichloromethane', 'ClCCl'),
    ('ethane-1,2-diol', 'OCCO'),
    ('ethanol', 'CCO'),
    ('ethoxybenzene', 'CCOc1ccccc1'),
    ('ethoxyethane', 'CCOCC'),
    ('ethyl acetate', 'CCOC(C)=O'),
    ('ethyl benzoate', 'CCOC(=O)C1=CC=CC=C1'),
    ('formamide', 'NC=O'),
    ('heavy water', '[2H]O[2H]'),
    ('heptan-1-ol', 'CCCCCCCO'),
    ('heptane', 'CCCCCCC'),
    ('hexamethylphosphoramide', 'CN(C)P(=O)(N(C)C)N(C)C'),
    ('hexan-1-ol', 'CCCCCCO'),
    ('hexane', 'CCCCCC'),
    ('methanedithione', 'S=C=S'),
    ('methanesulfinylmethane', 'CS(C)=O'),
    ('methanol', 'CO'),
    ('methoxybenzene', 'COc1ccccc1'),
    ('n,n-dimethylacetamide', 'CN(C)C(C)=O'),
    ('n,n-dimethylformamide', 'CN(C)C=O'),
    ('nitrobenzene', '[O-][N+](=O)C1=CC=CC=C1'),
    ('nitromethane', 'C[N+]([O-])=O'),
    ('octan-1-ol', 'CCCCCCCCO'),
    ('oxolan-2-one', 'O=C1CCCO1'),
    ('oxolane', 'C1CCOC1'),
    ('p-xylene', 'Cc1ccc(C)cc1'),
    ('pentan-1-ol', 'CCCCCO'),
    ('phenylmethanol', 'OCC1=CC=CC=C1'),
    ('piperidine', 'C1CCNCC1'),
    ('propan-1-ol', 'CCCO'),
    ('propan-2-ol', 'CC(C)O'),
    ('propan-2-one', 'CC(C)=O'),
    ('propanenitrile', 'CCC##N'),
    ('pyridine', 'C1=CC=NC=C1'),
    ('tetrachloromethane', 'ClC(Cl)(Cl)Cl'),
    ('tetrahydrothiophene 1,1-dioxide', 'O=S1(=O)CCCC1'),
    ('toluene', 'CC1=CC=CC=C1'),
    ('trichloromethane', 'ClC(Cl)Cl'),
    ('water', 'O'),
    # DON'T TOUCH ORDER. ONLY APPEND.
    ))


__all__ = ['Conditions', 'DictToConditions', 'ConditionsToDataFrame']
