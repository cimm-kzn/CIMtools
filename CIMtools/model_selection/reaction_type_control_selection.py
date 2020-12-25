# -*- coding: utf-8 -*-
#
#  Copyright 2019 Assima Rakhimbekova <asima.astana@outlook.com>
#  Copyright 2019, 2020 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from CGRtools.containers import ReactionContainer
from numpy import hstack
from sklearn.model_selection import KFold
from sklearn.utils import _safe_indexing, indexable
from ..applicability_domain import ReactionTypeControl
from ..metrics import balanced_accuracy_score_with_ad, rmse_score_with_ad
from ..utils import iter2array


def rtc_env_selection(X, y, data, envs, reg_model, score):
    """
    Function for finding the best number of neighbours in ReactionTypeControl method.

    All AD’s model hyperparameters were selected based on internal cross-validation using training set.
    The hyperparameters of the AD definition approach have been optimized in the cross-validation,
    where metrics RMSE_AD or BA_AD were used as maximized scoring functions.

    :param X: array-like or sparse matrix, shape (n_samples, n_features)
            The input samples. Internally, it will be converted to
            ``dtype=np.float32`` and if a sparse matrix is provided
             to a sparse ``csr_matrix``.
    :param y: array-like, shape = [n_samples] or [n_samples, n_outputs]
             The target values (real numbers in regression).
    :param data: after read rdf file
    :param envs: list or tuple. Numbers of neighbours.
    :param reg_model: estimator
    :param score: 'ba_ad' or 'rmse_ad'
    :return: int
    """
    data = iter2array(data, dtype=ReactionContainer)
    X, y, data = indexable(X, y, data)

    if not isinstance(envs, (list, tuple)):
        raise ValueError('envs must be list or tuple.')
    if reg_model is None:
        raise ValueError('Model is not defined.')
    if score not in ('ba_ad', 'rmse_ad'):
        raise ValueError('Invalid value for score. Allowed string values are "ba_ad", "rmse_ad".')

    cv = KFold(n_splits=5, shuffle=True, random_state=1)
    score_value = 0
    env_value = 0
    for env in envs:
        Y_pred, Y_true, AD = [], [], []
        for train_index, test_index in cv.split(X):
            x_train = _safe_indexing(X, train_index)
            x_test = _safe_indexing(X, test_index)
            y_train = _safe_indexing(y, train_index)
            y_test = _safe_indexing(y, test_index)
            data_train = _safe_indexing(data, train_index)
            data_test = _safe_indexing(data, test_index)
            Y_pred.append(reg_model.fit(x_train, y_train).predict(x_test))
            Y_true.append(y_test)
            AD.append(ReactionTypeControl(env=env).fit(data_train).predict(data_test))
        if score == 'ba_ad':
            val = balanced_accuracy_score_with_ad(Y_true=hstack(Y_true), Y_pred=hstack(Y_pred), AD=hstack(AD))
        else:
            val = rmse_score_with_ad(Y_true=hstack(Y_true), Y_pred=hstack(Y_pred), AD=hstack(AD))
        if val >= score_value:
            score_value = val
            env_value = env
    return env_value


__all__ = ['rtc_env_selection']
