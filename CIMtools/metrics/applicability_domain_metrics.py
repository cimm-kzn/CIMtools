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
from numpy import sqrt, hstack, unique
from sklearn.metrics import balanced_accuracy_score, mean_squared_error
from sklearn.model_selection import KFold
from sklearn.utils import safe_indexing
from CIMtools.applicability_domain import ReactionTypeControl

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


def optimal_env(X, y, data, envs, reg_model, score):
    cv = KFold(n_splits=5, shuffle=True, random_state=1)
    score_value = 0
    env_value = 0
    for env in envs:
        Y_pred, Y_true, AD = [], [], []
        for train_index, test_index in cv.split(X):
            x_train = safe_indexing(X, train_index)
            x_test = safe_indexing(X, test_index)
            y_train = safe_indexing(y, train_index)
            y_test = safe_indexing(y, test_index)
            data_train = safe_indexing(data, train_index)
            data_test = safe_indexing(data, test_index)
            Y_pred.append(reg_model.fit(x_train, y_train).predict(x_test))
            Y_true.append(y_test)
            AD.append(ReactionTypeControl(env=env).fit(data_train).predict(data_test))
        AD_stack = hstack(AD)
        AD_ = unique(AD_stack)
        for z in AD_:
            AD_new = AD_stack <= z
        if score == 'ba_ad':
            val = balanced_accuracy_score_with_ad(Y_true=hstack(Y_true), Y_pred=hstack(Y_pred), AD=AD_new)
        elif score == 'rmse_ad':
            val = rmse_score_with_ad(Y_true=hstack(Y_true), Y_pred=hstack(Y_pred), AD=AD_new)
        if val >= score_value:
            score_value = val
            env_value = env
    return env_value


__all__ = ['balanced_accuracy_score_with_ad', 'rmse_score_with_ad', 'optimal_env']
