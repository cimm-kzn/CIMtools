# -*- coding: utf-8 -*-
#
#  Copyright 2018, 2019 Ramil Nugmanov <stsouko@live.ru>
#  Copyright 2019 Ravil Mukhametgaleev <sonic-mc@mail.ru>
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
from functools import partial
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
            k = hyphens_replace(k)
            k = single_quotes_replace(k)
            k = k.replace("''", '"')  # double quotes from single
            try:
                k = known_solvents[k.lower()]
            except KeyError:
                raise KeyError(f'unknown solvent: {k}')
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


known_solvents = (
    ('1-phenylethan-1-one', 'CC(=O)c1ccccc1', 'acetophenone', 'methyl phenyl ketone', 'phenylethanone'),
    ('1,2-dichloroethane', 'ClCCCl', 'ethylene dichloride', '1,2-dca', 'dce', 'ethane dichloride', 'dutch liquid',
     'dutch oil', 'freon 150', 'freon-150'),
    ('1,2-dimethoxyethane', 'COCCOC', 'ethane-1,2-diyl dimethyl ether', 'dme', 'glyme',
     'ethylene glycol dimethyl ether', 'monoglyme', 'dimethyl glycol', 'dimethyl cellosolve'),
    ('propane-1,2,3-triol', 'OCC(O)CO', '1,2,3-propanetriol', 'glycerin', 'glycerine', 'propanetriol',
     '1,2,3-trihydroxypropane'),
    ('1,3-dimethylbenzene', 'CC1=CC(C)=CC=C1', '1,3-xylene', '3-xylene', 'isoxylene', 'm-xylene', 'm-xylol'),
    ('1,3,5-trimethylbenzene', 'CC1=CC(C)=CC(C)=C1', 'mesitylene', 'sym-trimethylbenzene'),
    ('1,4-dimethylbenzene', 'Cc1ccc(C)cc1', 'p-xylene', '1,4-xylene', 'p-dimethylbenzene', 'p-xylol', 'p-methyltoluene',
     'paraxylene', 'chromar', 'scintillar', '4-methyltoluene', 'nsc 72419'),
    ('1,4-dioxane', 'C1COCCO1', '1,4-dioxacyclohexane', '[1,4]dioxane', 'p-dioxane', '[6]-crown-2',
     'diethylene dioxide', 'diethylene ether', 'dioxan'),
    ('1,4-epoxybutane', 'C1CCOC1', 'oxolane', 'oxacyclopentane', 'tetrahydrofuran', 'butylene oxide',
     'cyclotetramethylene oxide', 'diethylene oxide', 'tetra-methylene oxide'),
    ('2-methylbutan-2-ol', 'CCC(C)(C)O', '2-methyl-2-butanol', 'tert-amyl alcohol', 't-amylol', 'taa',
     'tert-pentyl alcohol', '2-methyl-2-butyl alcohol', 't-pentylol', 'amylene hydrate', 'dimethylethylcarbinol'),
    ('2-methylpropan-1-ol', 'CC(C)CO', 'isobutyl alcohol', 'iba', '2-methyl-1-propanol', '2-methylpropyl alcohol',
     'isopropylcarbinol'),
    ('2-methylpropan-2-ol', 'CC(C)(C)O', 't-butyl alcohol', 'tert-butanol', 't-butanol', 'trimethylcarbinol',
     '2-methyl-2-propanol', '2m2p'),
    ('2,2,4-trimethylpentane', 'CC(C)CC(C)(C)C', 'isobutyltrimethylmethane', 'isooctane'),
    ('3-methylbutan-1-ol', 'CC(C)CCO', '3-methyl-1-butanol', '3-methylbutanol', 'i-amyl alcohol', 'isoamyl alcohol',
     'isobutyl carbinol', 'isopentanol', 'isopentyl alcohol'),
    ('acetic acid', 'CC(O)=O', 'ethanoic acid', 'glacial acetic acid'),
    ('acetonitrile', 'CC##N', 'ethanenitrile', 'cyanomethane', 'ethyl nitrile', 'methanecarbonitrile',
     'methyl cyanide'),
    ('benzene', 'C1=CC=CC=C1', 'benzol', '[6]annulene'),
    ('benzenecarbonitrile', 'N#Cc1ccccc1', 'benzonitrile', 'cyanobenzene', 'phenyl cyanide'),
    ('bromobenzene', 'BrC1=CC=CC=C1', 'phenyl bromide', 'bromobenzol', 'monobromobenzene'),
    ('butan-1-ol', 'CCCCO', '1-butanol', '1-butyl alcohol', 'butanol', 'butyl alcohol', 'n-butanol',
     'n-butyl alcohol', 'butalcohol', 'butyl hydrate', 'butylic alcohol', 'butyralcohol', 'butyric alcohol',
     'butyryl alcohol', '1-hydroxybutane', 'n-propylcarbinol'),
    ('butan-2-ol', 'CCC(C)O', '2-butanol', 's-butanol', 'sec-butanol', 'sec-butyl alcohol', '2-butyl alcohol'),
    ('butan-2-one', 'CCC(C)=O', '2-butanone', 'butanone', 'ethyl methyl ketone', 'mek', 'methyl ethyl ketone'),
    ('chlorobenzene', 'ClC1=CC=CC=C1', 'benzene chloride', 'monochlorobenzene', 'phenyl chloride', 'chlorobenzol',
     'mcb'),
    ('cyclohexane', 'C1CCCCC1', 'hexahydrobenzene', 'hexamethylene'),
    ('deuterium oxide', '[2H]O[2H]', 'heavy water', 'dideuterium monoxide', 'water-d2'),
    ('dichloromethane', 'ClCCl',  'methylene chloride', 'methylene dichloride', 'methylenechloride'),
    ('ethane-1,2-diol', 'OCCO', 'ethylene glycol', '1,2-ethanediol', 'ethylene alcohol', 'hypodicarbonous acid',
     'monoethylene glycol', '1,2-dihydroxyethane'),
    ('ethanol', 'CCO', 'absolute alcohol', 'alcohol', 'cologne spirit', 'drinking alcohol', 'ethylic alcohol', 'etoh',
     'ethyl alcohol', 'ethyl hydrate', 'ethyl hydroxide', 'ethylol', 'grain alcohol', 'hydroxyethane',
     'methylcarbinol'),
    ('ethoxybenzene', 'CCOc1ccccc1', 'ethyl phenyl ether', 'phenetole', 'phenyl ethyl ether'),
    ('ethoxyethane', 'CCOCC', '1,1-oxybisethane', '3-oxapentane', 'diethyl ether', 'ether', 'diethyl ether', 'dether',
     'ethyl ether', 'ethyl oxide', '3-oxapentane', 'diethyl oxide', 'solvent ether', 'sulfuric ether'),
    ('ethyl acetate', 'CCOC(C)=O', 'ethyl ethanoate', 'acetic ester', 'acetic ether', 'ethyl ester of acetic acid'),
    ('ethyl benzoate', 'CCOC(=O)C1=CC=CC=C1', 'ethyl benzenecarboxylate', 'ethyl phenylformate'),
    ('formamide', 'NC=O', 'methanamide', 'carbamaldehyde', 'formyl amide', 'methanamide'),
    ('heptan-1-ol', 'CCCCCCCO', '1-heptanol', 'enenthic alcohol', 'heptanol', 'heptyl alcohol', 'n-heptyl alcohol'),
    ('heptane', 'CCCCCCC', 'n-heptane', 'septane'),
    ('hexamethylphosphoramide', 'CN(C)P(=O)(N(C)C)N(C)C', 'hexamethylphosphoric acid triamide',
     'hexamethylphosphoric triamide', 'hexamethylphosphorous triamide', 'hmpa',
     '''N,N,N',N',N",N"-hexamethylphosphoric triamide''', 'N-bis(dimethylamino)phosphoryl-N-methyl-methanamine',
     'tris(dimethylamino)phosphine oxide'),
    ('hexan-1-ol', 'CCCCCCO', '1-hexanol', 'hexanol', 'hexanol-1', 'hexyl alcohol', 'n-hexanol', 'n-hexyl alcohol'),
    ('hexane', 'CCCCCC', 'n-hexane', 'sextane'),
    ('methanedithione', 'S=C=S', 'carbon disulfide', 'carbon disulphide'),
    ('methanesulfinylmethane', 'CS(C)=O', 'dimethyl sulfoxide', 'dimethyl(oxido)sulfur', 'methylsulfinylmethane',
     'methyl sulfoxide', '(methanesulfinyl)methane'),
    ('methanol', 'CO', 'carbinol', 'columbian spirits', 'hydroxymethane', 'methyl alcohol', 'methyl hydrate',
     'methyl hydroxide', 'methylic alcohol', 'methylol', 'pyroligneous spirit', 'wood alcohol', 'wood naphtha',
     'wood spirit'),
    ('methoxybenzene', 'COc1ccccc1', 'anisol', 'anisole', 'methyl phenyl ether', 'phenoxymethane'),
    ('N,N-dimethylacetamide', 'CN(C)C(C)=O', 'acetyl dimethylamine', 'dimethylacetamide', 'dma',
     'N,N-dimethylethanamide'),
    ('N,N-dimethylformamide', 'CN(C)C=O', 'dimethylaminoformaldehyde', 'dimethylformamide', 'dmf',
     'N,N-dimethylmethanamide'),
    ('nitrobenzene', '[O-][N+](=O)C1=CC=CC=C1', 'nitrobenzol', 'oil of mirbane'),
    ('nitromethane', 'C[N+]([O-])=O', 'nitrocarbol'),
    ('octan-1-ol', 'CCCCCCCCO', '1-octanol', 'n-octanol', 'capryl alcohol', 'octyl alcohol'),
    ('oxolan-2-one', 'O=C1CCCO1', '1,4-butanolide', '2(3h)-dihydrofuranone', '4-butyrolactone',
     '4-hydroxybutanoic acid lactone', '4-hydroxybutyric acid g-lactone', 'butyrolactone', 'dihydro-2(3h)-furanone',
     'γ-butyrolactone', 'gamma-butyrolactone', 'g-butyrolactone', 'gamma butyrolactone'),
    ('pentan-1-ol', 'CCCCCO', '1-pentanol', 'amyl alcohol', 'butyl carbinol', 'n-amyl alcohol', 'n-pentanol',
     'n-pentyl alcohol', 'pentanol'),
    ('phenylmethanol', 'OCC1=CC=CC=C1', 'phenylcarbinol', 'benzenemethanol', 'benzyl alcohol'),
    ('piperidine', 'C1CCNCC1', 'azinane', 'hexahydropyridine', 'pentamethyleneamine', 'azacyclohexane', 'azinane'),
    ('propan-1-ol', 'CCCO', '1-propanol', 'n-propanol', 'n-propyl alcohol', 'propanol', 'propyl alcohol', 'n-proh',
     'ethylcarbinol', '1-hydroxypropane', 'propionic alcohol', 'propionyl alcohol', 'propionylol', 'propyl alcohol',
     'propylic alcohol', 'propylol'),
    ('propan-2-ol', 'CC(C)O', '2-propanol', 'isopropanol', 'isopropyl alcohol', 'rubbing alcohol', 'sec-propyl alcohol',
     's-propanol', 'iproh', 'i-proh', 'dimethyl carbinol', 'ipa'),
    ('propan-2-one', 'CC(C)=O', 'acetone', 'dimethyl ketone', 'dimethyl carbonyl', 'β-ketopropane', 'propanone',
     '2-propanone', 'dimethyl formaldehyde', 'pyroacetic spirit (archaic)', 'Ketone propane'),
    ('propanenitrile', 'CCC##N', 'cyanoethane', 'ethyl cyanide', 'propionitrile', 'propylnitrile', 'propiononitrile'),
    ('pyridine', 'C1=CC=NC=C1', 'azinine', 'azine', 'azabenzene'),
    ('tetrachloromethane', 'ClC(Cl)(Cl)Cl', 'carbon tetrachloride', 'benziform', 'benzinoform', 'carbon chloride',
     'carbon tet', 'freon-10', 'freon 10', 'refrigerant-10', 'halon-104', 'methane tetrachloride',
     'methyl tetrachloride', 'perchloromethane', 'tetraform', 'tetrasol'),
    ('tetrahydrothiophene 1,1-dioxide', 'O=S1(=O)CCCC1', 'sulfolane', 'tetrahydrothiophene dioxide',
     'tetramethylene sulfone', 'tetramethylene sulphone', 'thiacyclopentane dioxide', 'thiolane 1,1-dioxide',
     'thiolane-1,1-dioxide'),
    ('toluene', 'CC1=CC=CC=C1', 'methylbenzene', 'toluol', 'methyl benzene', 'phenyl methane', 'anisen'),
    ('trichloromethane', 'ClC(Cl)Cl', 'chloroform', 'methane trichloride', 'methyl trichloride', 'methenyl trichloride',
     'tcm', 'freon 20', 'freon-20', 'refrigerant-20', 'r-20', 'un 1888'),
    ('water', 'O')
)


def multi_replace(string: str, patterns, replacement: str):
    for p in patterns:
        string = string.replace(p, replacement)
    return string


single_quotes_replace = partial(multi_replace, patterns=('‘', '’'), replacement="'")
hyphens_replace = partial(multi_replace, patterns=('‐', '‑', '‒', '–', '—', '―', '₋', '−'), replacement='-')


tmp = {}
solvent_smiles = {}
for n, s, *o in known_solvents:
    solvent_smiles[n] = s
    low_n = n.lower()
    tmp[low_n] = n

    for x in o:
        tmp[x.lower()] = n

known_solvents = tmp


__all__ = ['Conditions', 'DictToConditions', 'ConditionsToDataFrame', 'known_solvents', 'solvent_smiles']
