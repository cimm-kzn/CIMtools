# -*- coding: utf-8 -*-
#
#  Copyright 2020 Assima Rakhimbekova <asima.astana@outlook.com>
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
from sklearn.base import BaseEstimator, clone
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.model_selection import KFold
from sklearn.preprocessing import StandardScaler
from sklearn.utils import safe_indexing
from sklearn.utils.validation import check_array, check_is_fitted
from ..metrics.applicability_domain_metrics import balanced_accuracy_score_with_ad, rmse_score_with_ad


class GPR_AD(BaseEstimator):
    """
    Gaussian Process Regression (GPR) assumes that the joint distribution of a real-valued property of
    chemical reactions and their descriptors is multivariate normal (Gaussian) with the elements of its covariance
    matrix computed by means of special covariance functions (kernels). For every reaction,
    a GPR model produces using the Bayes’ theorem a posterior conditional distribution (so-called prediction density)
    of the reaction property given the vector of reaction descriptors.
    The prediction density has normal (Gaussian) distribution with the mean corresponding to predicted value
    of the property and the variance corresponding to prediction confidence [1].
    If the variance is greater than a predefined threshold σ*, the chemical reaction is considered as
    X-outlier (out of AD)
    """
    def __init__(self, threshold='cv', score='ba_ad', gpr_model=None):
        self.threshold = threshold
        self.score = score
        self.gpr_model = gpr_model
        if threshold != 'cv' and not isinstance(threshold, float):
            raise ValueError('Invalid value for threshold. Allowed string value is "cv".')
        if score not in ('ba_ad', 'rmse_ad'):
            raise ValueError('Invalid value for score. Allowed string values are "ba_ad", "rmse_ad".')

    def fit(self, X, y=None):
        """
        Model building and threshold searching
        During training, a model is built and a ariance threshold σ* is found by which the object is considered to
        belong to the applicability domain of the model.
        X : array-like or sparse matrix, shape (n_samples, n_features)
            The input samples. Use ``dtype=np.float32`` for maximum
            efficiency.
        y : array-like, shape = [n_samples] or [n_samples, n_outputs]
            The target values (real numbers in regression).

        Returns
        -------
        self : object
        """
        if y is None:
            raise ValueError("Y must be specified to find the optimal threshold.")
        y = check_array(y, accept_sparse='csc', ensure_2d=False, dtype=None)
        # Check that X have correct shape
        X = check_array(X)

        self.scaler = StandardScaler()
        y_gpr = self.scaler.fit_transform(y.reshape(-1, 1))

        if self.gpr_model is None:
            gpr_model = GaussianProcessRegressor(random_state=1)
        else:
            gpr_model = clone(self.gpr_model)

        self.AD_gpr = gpr_model.fit(X, y_gpr)

        if self.threshold == 'cv':
            cv = KFold(n_splits=5, random_state=1, shuffle=True)
            gpr_model_int = clone(gpr_model)

            self.threshold_value = 0
            score_value = 0
            Y_pred, Y_true, AD = [], [], []
            for train_index, test_index in cv.split(X):
                x_train = safe_indexing(X, train_index)
                x_test = safe_indexing(X, test_index)
                y_train = safe_indexing(y_gpr, train_index)
                y_test = safe_indexing(y_gpr, test_index)
                gpr_model_int.fit(x_train, y_train)
                y_pred, y_var = gpr_model_int.predict(x_test, return_std=True)

                Y_pred.append(y_pred.flatten())
                Y_true.append(y_test.flatten())
                AD.append(y_var)
            AD_stack = hstack(AD)
            AD_ = unique(AD_stack)
            for z in AD_:
                AD_new = AD_stack <= z
                if self.score == 'ba_ad':
                    val = balanced_accuracy_score_with_ad(Y_true=hstack(Y_true), Y_pred=hstack(Y_pred), AD=AD_new)
                elif self.score == 'rmse_ad':
                    val = rmse_score_with_ad(Y_true=hstack(Y_true), Y_pred=hstack(Y_pred), AD=AD_new)
                if val >= score_value:
                    score_value = val
                    self.threshold_value = z
        else:
            self.threshold_value = self.threshold
        return self

    def predict(self, X):
        """Predict inside or outside AD for X.

        Parameters
        ----------
        X : array-like or sparse matrix, shape (n_samples, n_features)
            The input samples. Internally, it will be converted to
            ``dtype=np.float32`` and if a sparse matrix is provided
            to a sparse ``csr_matrix``.

        Returns
        -------
        ad : array of shape = [n_samples]
            Array contains True (reaction in AD) and False (reaction residing outside AD).
        """
        # Check is fit had been called
        check_is_fitted(self, ['AD_gpr', 'threshold_value'])
        # Check that X have correct shape
        X = check_array(X)
        Y_pred_GPR, Y_var = self.AD_gpr.predict(X, return_std=True)

        return Y_var <= self.threshold_value


__all__ = ['GPR_AD']
