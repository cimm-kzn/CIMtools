# -*- coding: utf-8 -*-
#
#  Copyright 2018-2020 Ramil Nugmanov <nougmanoff@protonmail.com>
#  Copyright 2021 Ilnura Kashafutdinova <ilnura12@mail.ru >
#  Copyright 2021 Assima Rakhimbekova <assima.rakhimbekova@outlook.com>
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
from itertools import compress
from pandas import DataFrame
from ..base import CIMtoolsTransformerMixin
from ..exceptions import ConfigurationError


class SolventVectorizer(CIMtoolsTransformerMixin):
    def __init__(self, polarizability_form1=True, polarizability_form2=True, permettivity_form1=True,
                 permettivity_form2=True, permettivity_form3=True, permettivity_form4=True,
                 permettivity_polarizability=True, alpha_kamlet_taft=True, beta_kamlet_taft=True, pi_kamlet_taft=True,
                 spp_katalan=True, sb_katalan=True, sa_katalan=True, boiling_point=False, melting_point=False):
        self.polarizability_form1 = polarizability_form1
        self.polarizability_form2 = polarizability_form2
        self.permettivity_form1 = permettivity_form1
        self.permettivity_form2 = permettivity_form2
        self.permettivity_form3 = permettivity_form3
        self.permettivity_form4 = permettivity_form4
        self.permettivity_polarizability = permettivity_polarizability
        self.alpha_kamlet_taft = alpha_kamlet_taft
        self.beta_kamlet_taft = beta_kamlet_taft
        self.pi_kamlet_taft = pi_kamlet_taft
        self.spp_katalan = spp_katalan
        self.sb_katalan = sb_katalan
        self.sa_katalan = sa_katalan
        self.boiling_point = boiling_point
        self.melting_point = melting_point
        self.__prepare_header()

    def __getstate__(self):
        return {k: v for k, v in super().__getstate__().items() if not k.startswith('_SolventVectorizer__')}

    def __setstate__(self, state):
        if not state:  # < 4.0.5 reverse compatibility
            state = {'polarizability_form1': True, 'polarizability_form2': True, 'permettivity_form1': True,
                     'permettivity_form2': True, 'permettivity_form3': True, 'permettivity_form4': True,
                     'permettivity_polarizability': True, 'alpha_kamlet_taft': True, 'beta_kamlet_taft': True,
                     'pi_kamlet_taft': True, 'spp_katalan': True, 'sb_katalan': True, 'sa_katalan': True}
        super().__setstate__(state)
        self.__prepare_header()

    def set_params(self, **params):
        super().set_params(**params)
        self.__prepare_header()

    def __prepare_header(self):
        header = []
        index = []
        if self.polarizability_form1:
            header.append('Polarizability (form1)')
            index.append(True)
        else:
            index.append(False)
        if self.polarizability_form2:
            header.append('Polarizability (form2)')
            index.append(True)
        else:
            index.append(False)
        if self.permettivity_form1:
            header.append('Permettivity (form1)')
            index.append(True)
        else:
            index.append(False)
        if self.permettivity_form2:
            header.append('Permettivity (form2)')
            index.append(True)
        else:
            index.append(False)
        if self.permettivity_form3:
            header.append('Permettivity (form3)')
            index.append(True)
        else:
            index.append(False)
        if self.permettivity_form4:
            header.append('Permettivity (form4)')
            index.append(True)
        else:
            index.append(False)
        if self.permettivity_polarizability:
            header.append('Permettivity-Polarizability')
            index.append(True)
        else:
            index.append(False)
        if self.alpha_kamlet_taft:
            header.append('alpha Kamlet-Taft')
            index.append(True)
        else:
            index.append(False)
        if self.beta_kamlet_taft:
            header.append('beta Kamlet-Taft')
            index.append(True)
        else:
            index.append(False)
        if self.pi_kamlet_taft:
            header.append('pi Kamlet-Taft')
            index.append(True)
        else:
            index.append(False)
        if self.spp_katalan:
            header.append('SPP Katalan')
            index.append(True)
        else:
            index.append(False)
        if self.sb_katalan:
            header.append('SB Katalan')
            index.append(True)
        else:
            index.append(False)
        if self.sa_katalan:
            header.append('SA Katalan')
            index.append(True)
        else:
            index.append(False)

        if not header:
            raise ConfigurationError('required at least one parameter')
        self.__header = header
        self.__index = index

    def get_feature_names(self):
        """Get feature names.

        Returns
        -------
        feature_names : list of strings
            Names of the features produced by transform.
        """
        return self.__header

    def transform(self, x):
        x = super().transform(x)
        return DataFrame([compress(described_solvents[x], self.__index) for x in x], columns=self.__header)

    _dtype = str


described_solvents = dict((
        ('1-phenylethan-1-one', (.31, .237, .845, .943, .458, .891, .109, .04, .49, .9, .9, .37, .04, 202, 20)),
        ('1,2-dichloroethane', (.266, .21, .757, .903, .431, .824, .091, .0, .1, .81, .89, .13, .03, 83.5, -35.5)),
        ('1,2-dimethoxyethane', (.23, .19, .67, .86, .4, .76, .08, .0, .41, .53, .79, .64, .0, 85, -58)),
        ('propane-1,2,3-triol', (.28, .22, .94, .98, .48, .96, .11, 1.21, .51, .62, .95, .31, .65, 290, 18.2)),
        ('1,3-dimethylbenzene', (.293, .226, .318, .583, .241, .412, .055, .0, .11, .47, .62, .16, .0, 139.1, -47.8)),
        ('1,3,5‐trimethylbenzene', (.294, .227, .318, .583, .241, .412, .055, .0, .13, .41, .58, .19, .0, 164.7, -44.7)),
        ('1,4-dioxane', (.253, .202, .287, .548, .223, .377, .045, .0, .37, .49, .7, .44, .0, 101.5, 11.8)),
        ('2-methylbutan-2-ol', (.245, .197, .614, .827, .381, .705, .075, .28, .93, .4, .83, .94, .1, 102.4, -9.1)),
        ('2-methylpropan-1-ol', (.24, .194, .84, .94, .456, .887, .088, .79, .84, .4, .83, .83, .31, 107.8, -108)),
        ('2-methylpropan-2-ol', (.234, .19, .793, .92, .442, .852, .084, .42, .93, .41, .83, .93, .15, 82.4, 25.4)),
        ('2,2,4-trimethylpentane', (.237, .191, .242, .49, .195, .324, .037, .0, .0, .04, .53, .04, .0, 99.2, -107.3)),
        ('3-methylbutan-1-ol', (.245, .197, .825, .934, .452, .876, .089, .84, .86, .4, .81, .86, .32, 131.1, -117.2)),
        ('acetic acid', (.226, .184, .632, .837, .387, .72, .071, 1.12, .45, .64, .78, .39, .69, 117.9, 16.6)),
        ('acetonitrile', (.21, .174, .921, .972, .479, .946, .083, .19, .4, .66, .9, .29, .04, 81.6, -43.8)),
        ('benzene', (.293, .227, .297, .559, .229, .388, .052, .0, .1, .55, .67, .12, .0, 80, 5.5)),
        ('benzenecarbonitrile', (.31, .24, .89, .96, .47, .93, .11, .0, .37, .9, .96, .28, .05, 191.1, -12.7)),
        ('bromobenzene', (.323, .244, .595, .815, .373, .688, .091, .0, .06, .79, .82, .19, .0, 156, -30.6)),
        ('butan-1-ol', (.241, .194, .846, .943, .458, .892, .089, .84, .84, .47, .84, .81, .34, 117.85, -89.8)),
        ('butan-2-ol', (.24, .193, .838, .94, .456, .886, .088, .69, .8, .4, .84, .89, .22, 99.5, -114.7)),
        ('butan-2-one', (.23, .187, .851, .945, .46, .895, .086, .06, .48, .6, .88, .52, .0, 79.5, -86.6)),
        ('chlorobenzene', (.304, .233, .606, .822, .377, .698, .088, .0, .07, .68, .82, .18, .0, -131.7, -45.2)),
        ('cyclohexane', (.255, .203, .254, .505, .202, .338, .041, .0, .0, .0, .56, .07, .0, 80.7, 6.6)),
        ('dichloromethane', (.254, .202, .726, .888, .42, .799, .085, .13, .1, .82, .88, .18, .04, 40, -95.1)),
        ('ethane-1,2-diol', (.259, .205, .924, .973, .48, .948, .099, .9, .52, .92, .93, .53, .72, 197.3, -13)),
        ('ethanol', (.22, .181, .887, .959, .47, .922, .085, .86, .75, .54, .85, .66, .4, 78.2, -114.1)),
        ('ethoxybenzene', (.3, .23, .52, .76, .34, .62, .08, .0, .3, .69, .74, .3, .0, 169.8, -29.5)),
        ('ethoxyethane', (.215, .177, .516, .762, .34, .615, .06, .0, .47, .24, .69, .56, .0, 34.6, -116.3)),
        ('ethyl acetate', (.228, .19, .63, .83, .38, .72, .071, .0, .45, .55, .8, .54, .0, 77.1, -83.6)),
        ('ethyl benzoate', (.297, .229, .625, .833, .385, .714, .088, .0, .41, .74, .84, .42, .0, 212, -34)),
        ('formamide', (.267, .211, .973, .991, .493, .982, .104, .71, .48, .97, .83, .41, .55, 220, 2.55)),
        ('deuterium oxide', (.203, .169, .963, .987, .49, .975, .083, 1.17, .47, 1.09, .96, .44, 1.06, 101.42, 3.81)),
        ('heptan-1-ol', (.255, .203, .775, .912, .437, .838, .089, .79, .82, .4, .8, .91, .3, 176.4, -34)),
        ('heptane', (.236, .19, .235, .479, .19, .315, .036, .0, .0, -.08, .53, .08, .0, 98.5, -90.6)),
        ('hexamethylphosphoramide', (.277, .217, .906, .967, .475, .935, .103, .0, 1.05, .87, .93, .81, .0, 232.5, 7.2)),
        ('hexan-1-ol', (.251, .201, .804, .925, .446, .86, .089, .8, .84, .4, .81, .88, .32, 157.6, -44.6)),
        ('hexane', (.227, .185, .227, .468, .185, .306, .034, .0, .0, -.11, .52, .06, .0, 68.7, -95.3)),
        ('methanedithione', (.355, .262, .348, .615, .258, .444, .068, .0, .07, .61, .59, .1, .0, 46, -111.5)),
        ('methanesulfinylmethane', (.283, .22, .938, .978, .484, .958, .107, .0, .76, 1.0, 1.0, .65, .07, 189, 18.5)),
        ('methanol', (.202, .168, .913, .969, .477, .941, .08, .98, .66, .6, .86, .55, .61, 64.7, -97.6)),
        ('methoxybenzene', (.3, .23, .52, .77, .34, .62, .08, .0, .32, .73, .82, .3, .08, 153.7, -37.5)),
        ('N,N-dimethylacetamide', (.26, .21, .92, .97, .48, .95, .1, .0, .76, .88, .97, .61, .03, 165, -20)),
        ('N,N-dimethylformamide', (.257, .205, .923, .973, .48, .947, .098, .0, .69, .88, .95, .61, .03, 153, -60.4)),
        ('nitrobenzene', (.319, .242, .918, .971, .479, .944, .116, .0, .3, .86, .97, .24, .06, 210.8, 5.7)),
        ('nitromethane', (.231, .188, .921, .972, .479, .946, .09, .22, .06, .75, .91, .24, .08, 101.1, -28.5)),
        ('octan-1-ol', (.257, .204, .757, .903, .431, .824, .088, .77, .81, .4, .79, .92, .3, 195.1, -15.5)),
        ('oxolan-2-one', (.26, .207, .927, .974, .481, .95, .099, .0, .49, .85, .99, .4, .06, 204, -43.3)),
        ('1,4-epoxybutane', (.245, .197, .687, .868, .407, .767, .08, .0, .55, .55, .84, .59, .0, 65, -108.3)),
        ('1,4-dimethylbenzene', (.29, .226, .302, .565, .232, .394, .052, .0, .12, .43, .62, .16, .0, 138.3, 13.2)),
        ('pentan-1-ol', (.247, .198, .811, .928, .448, .866, .089, .84, .86, .4, .82, .86, .32, 137.9, -78.9)),
        ('phenylmethanol', (.313, .238, .796, .921, .443, .854, .106, .6, .52, .98, .89, .46, .41, 205.3, -15.2)),
        ('piperidine', (.27, .213, .62, .83, .383, .71, .081, .0, 1.04, .3, .74, .93, .0, 106, -7)),
        ('propan-1-ol', (.234, .189, .866, .951, .464, .907, .088, .84, .9, .52, .85, .78, .37, 97.2, -126.1)),
        ('propan-2-ol', (.229, .186, .863, .95, .463, .904, .086, .76, .84, .48, .85, .83, .28, 82.9, -89.5)),
        ('propan-2-one', (.218, .179, .867, .951, .464, .907, .083, .08, .48, .62, .88, .48, .0, 55.75, -95.55)),
        ('propanenitrile', (.222, .182, .901, .965, .474, .932, .086, .0, .37, .64, .88, .37, .03, 97.1, -92.8)),
        ('pyridine', (.298, .229, .799, .923, .444, .856, .102, .0, .64, .87, .92, .58, .03, 115.2, -41.6)),
        ('tetrachloromethane', (.272, .214, .292, .554, .226, .383, .048, .0, .1, .21, .63, .04, .0, 76.8, -23)),
        ('tetrahydrothiophene 1,1-dioxide', (.285, .222, .934, .977, .483, .955, .107, .0, .39, .9, 1.0, .37, .05, 285, 27.6)),
        ('toluene', (.291, .226, .315, .58, .24, .408, .054, .0, .11, .49, .66, .13, .0, 110.6, -94.9)),
        ('trichloromethane', (.265, .209, .565, .796, .361, .66, .076, .2, .1, .58, .79, .07, .05, 61.1, -63.6)),
        ('water', (.205, .17, .963, .987, .49, .975, .084, 1.17, .47, 1.09, .96, .03, 1.06, 100, 0))))

__all__ = ['SolventVectorizer']
