# -*- coding: utf-8 -*-
#
#  Copyright 2019 Assima Rakhimbekova <asima.astana@outlook.com>
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
from numpy import sqrt
from sklearn.metrics import confusion_matrix


def function(metric, Y_true, Y_pred, AD):
    if metric == 'ba':
        tn, fp, fn, tp = confusion_matrix(Y_pred, AD).ravel()
        sen_TPR = (float(tp))/(tp+fn)
        spc_TNR = (float(tn)/(tn+fp))
        val = ((sen_TPR+spc_TNR)/2)
    else:
        AD_out = AD - 1
        AD_out_n = (AD_out == -1)
        s_n = sum(AD)
        s_out_n = sum(AD_out_n)
        if s_n:
            RMSE_AD = sqrt((sum(map(lambda x: (((x[0] - x[1]) ** 2) * x[2]), zip(Y_pred, Y_true, AD)))) / s_n)
        else:
            RMSE_AD = 0
        if s_out_n:
            RMSE_AD_out_n = sqrt(
                (sum(map(lambda x: (((x[0] - x[1]) ** 2) * x[2]), zip(Y_pred, Y_true, AD_out_n)))) / s_out_n)
        else:
            RMSE_AD_out_n = 0
        val = RMSE_AD_out_n - RMSE_AD
    return val
