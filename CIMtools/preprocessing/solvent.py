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
from pandas import DataFrame
from sklearn.base import BaseEstimator
from ..base import CGRtoolsTransformerMixin
from ..conditions_container import Conditions


class DominantSolvent(BaseEstimator, CGRtoolsTransformerMixin):
    def __init__(self):
        pass  # need for sklearn

    def transform(self, x):
        return DataFrame([described_solvents[x.solvents[0][0]] for x in super().transform(x)], columns=header)

    _dtype = Conditions


header = ('polarizability_form1', 'polarizability_form2', 'permettivity_form1', 'permettivity_form2',
          'permettivity_form3', 'permettivity_form4', 'permettivity_polarizability', 'alpha_kamlet-taft',
          'beta_Kamlet-Taft', 'pi*_Kamlet-taft', 'SPP_Katalan', 'SB_Katalan', 'SA_Katalan')
described_solvents = dict((
    ('1,3-xylene', (.293, .226, .318, .583, .241, .412, .055, .0, .11, .47, .62, .16, .0)),
    ('anisole', (.303, .232, .524, .767, .344, .623, .08, .0, .32, .73, .82, .30, .08)),
    ('1,2-dichloroethane', (.266, .21, .757, .903, .431, .824, .091, .0, .10, .81, .89, .13, .03)),
    ('benzene', (.293, .227, .297, .559, .229, .388, .052, .0, .10, .55, .67, .12, .0)),
    ('formamide', (.267, .211, .973, .991, .493, .982, .104, .71, .48, .97, .83, .41, .55)),
    ('ethoxyethane', (.215, .177, .516, .762, .34, .615, .06, .0, .47, .24, .69, .56, .0)),
    ('ethyl benzoate', (.297, .229, .625, .833, .385, .714, .088, .0, .41, .74, .84, .42, .0)),
    ('heptan-1-ol', (.255, .203, .775, .912, .437, .838, .089, .79, .82, .40, .80, .91, .30)),
    ('3-methylbutan-1-ol', (.245, .197, .825, .934, .452, .876, .089, .84, .86, .40, .81, .86, .32)),
    ('bromobenzene', (.323, .244, .595, .815, .373, .688, .091, .0, .06, .79, .82, .19, .0)),
    ('butan-1-ol', (.241, .194, .846, .943, .458, .892, .089, .84, .84, .47, .84, .81, .34)),
    ('butan-2-ol', (.24, .193, .838, .94, .456, .886, .088, .69, .80, .40, .84, .89, .22)),
    ('pentan-1-ol', (.247, .198, .811, .928, .448, .866, .089, .84, .86, .40, .82, .86, .32)),
    ('propan-2-ol', (.229, .186, .863, .95, .463, .904, .086, .76, .84, .48, .85, .83, .28)),
    ('cyclohexane', (.255, .203, .254, .505, .202, .338, .041, .0, .0, .0, .56, .07, .0)),
    ('propan-2-one', (.218, .179, .867, .951, .464, .907, .083, .08, .48, .62, .88, .48, .0)),
    ('chlorobenzene', (.304, .233, .606, .822, .377, .698, .088, .0, .07, .68, .82, .18, .0)),
    ('oxolane', (.245, .197, .687, .868, .407, .767, .08, .0, .55, .55, .84, .59, .0)),
    ('methanedithione', (.355, .262, .348, .615, .258, .444, .068, .0, .07, .61, .59, .10, .0)),
    ('pyridine', (.298, .229, .799, .923, .444, .856, .102, .0, .64, .87, .92, .58, .03)),
    ('acetic acid', (.226, .184, .632, .837, .387, .72, .071, 1.12, .45, .64, .78, .39, .69)),
    ('butan-2-one', (.23, .187, .851, .945, .46, .895, .086, .06, .48, .60, .88, .52, .0)),
    ('acetonitrile', (.21, .174, .921, .972, .479, .946, .083, .19, .40, .66, .90, .29, .04)),
    ('2,2,4-trimethylpentane', (.237, .191, .242, .49, .195, .324, .037, .0, .0, .04, .53, .04, .0)),
    ('2-methylpropan-1-ol', (.24, .194, .84, .94, .456, .887, .088, .79, .84, .40, .83, .83, .31)),
    ('2-methylpropan-2-ol', (.234, .19, .793, .92, .442, .852, .084, .42, .93, .41, .83, .93, .15)),
    ('methanol', (.202, .168, .913, .969, .477, .941, .08, .98, .66, .60, .86, .55, .61)),
    ('propan-1-ol', (.234, .189, .866, .951, .464, .907, .088, .84, .90, .52, .85, .78, .37)),
    ('methanesulfinylmethane', (.283, .22, .938, .978, .484, .958, .107, .0, .76, 1.00, 1.00, .65, .07)),
    ('dichloromethane', (.254, .202, .726, .888, .42, .799, .085, .13, .10, .82, .88, .18, .04)),
    ('ethanol', (.22, .181, .887, .959, .47, .922, .085, .86, .75, .54, .85, .66, .40)),
    ('water', (.205, .17, .963, .987, .49, .975, .084, 1.17, .47, 1.09, .96, .03, 1.06)),
    ('2-methylbutan-2-ol', (.245, .197, .614, .827, .381, .705, .075, .28, .93, .40, .83, .94, .10)),
    ('butanoate', (.227, .185, .625, .833, .385, .714, .071, .0, .45, .55, .80, .54, .0)),  # ERROR?
    ('1-phenylethan-1-one', (.31, .237, .845, .943, .458, .891, .109, .04, .49, .90, .90, .37, .04)),
    ('1,3,5‐trimethylbenzene', (.294, .227, .318, .583, .241, .412, .055, .0, .13, .41, .58, .19, .0)),
    ('ethyl acetate', (.227, .185, .625, .833, .385, .714, .071, .0, .45, .55, .80, .54, .0)),  # ERROR?
    ('1,4-dioxane', (.253, .202, .287, .548, .223, .377, .045, .0, .37, .49, .70, .44, .0)),
    ('tetrachloromethane', (.272, .214, .292, .554, .226, .383, .048, .0, .10, .21, .63, .04, .0)),
    ('ethane-1,2-diol', (.259, .205, .924, .973, .48, .948, .099, .90, .52, .92, .93, .53, .72)),
    ('trichloromethane', (.265, .209, .565, .796, .361, .66, .076, .20, .10, .58, .79, .07, .05)),
    ('propanenitrile', (.222, .182, .901, .965, .474, .932, .086, .0, .37, .64, .88, .37, .03)),
    ('phenylmethanol', (.313, .238, .796, .921, .443, .854, .106, .60, .52, .98, .89, .46, .41)),
    ('oxolan-2-one', (.26, .207, .927, .974, .481, .95, .099, .0, .49, .85, .99, .40, .06)),
    ('hexan-1-ol', (.251, .201, .804, .925, .446, .86, .089, .80, .84, .40, .81, .88, .32)),
    ('hexane', (.227, .185, .227, .468, .185, .306, .034, .0, .0, -0.11, .52, .06, .0)),
    ('octan-1-ol', (.257, .204, .757, .903, .431, .824, .088, .77, .81, .40, .79, .92, .3)),
    ('nitromethane', (.231, .188, .921, .972, .479, .946, .09, .22, .06, .75, .91, .24, .08)),
    ('nitrobenzene', (.319, .242, .918, .971, .479, .944, .116, .0, .30, .86, .97, .24, .06)),
    ('toluene', (.291, .226, .315, .58, .24, .408, .054, .0, .11, .49, .66, .13, .0)),
    ('heavy water', (.203, .169, .963, .987, .49, .975, .083, 1.17, .47, 1.09, .96, .44, 1.06)),
    ('tetrahydrothiophene 1,1-dioxide', (.285, .222, .934, .977, .483, .955, .107, .0, .39, .90, 1., .37, .05)),
    ('hexamethylphosphoramide', (.277, .217, .906, .967, .475, .935, .103, .0, 1.05, .87, .93, .81, .0)),
    ('n,n-dimethylformamide', (.257, .205, .923, .973, .48, .947, .098, .0, .69, .88, .95, .61, .03)),
    ('n,n-dimethylacetamide', (.26, .21, .92, .97, .48, .95, .10, .0, .76, .88, .97, .61, .03)),
    ('1,2,3-propanetriol', (.28, .22, .94, .98, .48, .96, .11, 1.21, .51, .62, .95, .31, .65)),
    ('1,2-dimethoxyethane', (.23, .19, .67, .86, .40, .76, .08, .0, .41, .53, .79, .64, .0))))


__all__ = ['DominantSolvent']
