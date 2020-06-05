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
    """
    The metric is used to optimize AD thresholds and their (hyper)parameters.
    This metric shows how well AD definition detects Y-outliers. First, property prediction errors are estimated
    in cross-validation for all reactions in a dataset. The reactions for which the absolute prediction error
    is higher than 3×RMSE are identified as Y-outliers, while the rest are considered as Y-inliers.
    Y-Outliers (poorly predicted) that are predicted by AD definition as X-outliers (outside AD) are called
    true outliers (TO), while Y-inliers predicted by AD definition as X-inliers (within AD) are called
    true inliers (TI). False outliers (FO) are Y-inliers that are wrongly predicted by the AD definition
    as X-outliers, while false inliers (FI) are Y-outliers that are wrongly predicted by the AD definition
    as X-inliers. The quality of outliers/inliers determination can be assessed using an analogue
    of the balanced accuracy.
    Parameters
    ----------
    Y_true: array-like, shape = [n_samples]
            The target values (real numbers in regression).
    Y_pred: array-like, shape = [n_samples]
            The predicted values of Y_true.
    AD: array-like, shape = [n_samples]
        Array contains True (reaction in AD) and False (reaction residing outside AD).

    Returns
    -------
    balanced_accuracy : float
    """
    AD_true = abs(Y_true - Y_pred) <= 3 * sqrt(mean_squared_error(Y_true, Y_pred))
    return balanced_accuracy_score(AD_true, AD)


def rmse_score_with_ad(Y_true, Y_pred, AD):
    """
    The metric is used to optimize AD thresholds and their (hyper)parameters.
    This metric is the difference between RMSE of property prediction for reactions outside AD and within AD.
    The metric was first proposed by Sahigata et al [1]. Negative values indicate that the reactions detected
    X-outliers (outside AD) are predicted better than X-inliers (within AD), thus highlighting some possible
    drawbacks in the definition of interpolation space. Its positive values indicate a reliable partition
    for the reactions detected as inside and outside AD and higher predictive performance within
    AD as compared to outside it. If no reactions are left inside or outside AD, then OIR is considered equal to 0.

    Parameters
    ----------
    Y_true: array-like, shape = [n_samples]
            The target values (real numbers in regression).
    Y_pred: array-like, shape = [n_samples]
            The predicted values of Y_true.
    AD: array-like, shape = [n_samples]
        Array contains True (reaction in AD) and False (reaction residing outside AD).
    Returns
    -------
    difference between RMSE of property prediction for reactions outside AD and within AD : float

    References
    -----------
    .. [1] Sahigara F., Mansouri K., Ballabio D., Mauri A., Consonni V. Todeschini R. Comparison of Different
           Approaches to Define the Applicability Domain of QSAR Models.  Molecules, 2012, vol. 17, pp. 4791-4810.
           doi: 10.3390/molecules17054791.
    """
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


__all__ = ['balanced_accuracy_score_with_ad', 'rmse_score_with_ad']
