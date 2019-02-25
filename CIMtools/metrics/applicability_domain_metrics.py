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
from sklearn.metrics import balanced_accuracy_score, mean_squared_error


def balanced_accuracy_score_with_ad(Y_true, Y_pred, AD):
    AD_true = abs(Y_true - Y_pred) <= 3 * sqrt(mean_squared_error(Y_true, Y_pred))
    return balanced_accuracy_score(AD_true, AD)

def rmse_score_with_ad(Y_true, Y_pred, AD):
    AD_out_n = ~AD
    s_n = AD.sum()
    s_out_n = AD_out_n.sum()
    if s_n:
        RMSE_AD = sqrt((sum((x - y) ** 2 for x, y, z in zip(Y_pred, Y_true, AD) if z)) / s_n)
    else:
        RMSE_AD = 0
    if s_out_n:
        RMSE_AD_out_n = sqrt((sum((x - y) ** 2 for x, y, z in zip(Y_pred, Y_true, AD_out_n) if z)) / s_out_n)
    else:
        RMSE_AD_out_n = 0
    return RMSE_AD_out_n - RMSE_AD
