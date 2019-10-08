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
from numpy import array, column_stack, eye, hstack, linalg, ones, unique
from sklearn.base import BaseEstimator, clone
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import KFold
from sklearn.utils import safe_indexing
from sklearn.utils.validation import check_array, check_is_fitted
from ..metrics.applicability_domain_metrics import balanced_accuracy_score_with_ad, rmse_score_with_ad


class Leverage(BaseEstimator):
    """ Distance-based method
    The model space can be represented by a two-dimensional matrix comprising n chemicals (rows) and
    k variables (columns), called the descriptor matrix (X). The leverage of a chemical provides a measure of the
    distance of the chemical from the centroid of X. Chemicals close to the centroid are less influential in model
    building than are extreme points. The leverages of all chemicals in the data set are generated by manipulating X
    according to Equation 1, to give the so-called Influence Matrix or Hat Matrix (H).

    H = X(XTX)–1 XT (Equation 1)

    where X is the descriptor matrix, XT is the transpose of X, and (A)–1 is the inverse of matrix A, where A = (XTX).

    The leverages or hat values (hi) of the chemicals (i) in the descriptor space are the diagonal elements of H,
    and can be computed by Equation 2.

    hii = xiT(XTX)–1 xi (Equation 2)

    where xi is the descriptor row-vector of the query chemical.
    A “warning leverage” (h*) is generally (!) fixed at 3p/n, where n is the number of training chemicals, and p
    the number of model variables plus one.

    A “warning leverage” can be found on internal cross-validation.

    A chemical with high leverage in the training set greatly influences
    the regression line: the fitted regression line is forced near to the observed value and its residual
    (observed-predicted value) is small, so the chemical does not appear to be an outlier, even though it may actually
    be outside the AD. In contrast, if a chemical in the test set has a hat value greater than the warning leverage h*,
    this means that the prediction is the result of substantial extrapolation and therefore may not be reliable.
    """
    def __init__(self, threshold='auto', score='ba_ad', reg_model=None):
        self.threshold = threshold
        self.score = score
        self.reg_model = reg_model
        if threshold not in ('auto', 'cv') and not isinstance(threshold, float):
            raise ValueError('Invalid value for threshold. Allowed string values are "auto", "cv".')
        if score not in ('ba_ad', 'rmse_ad'):
            raise ValueError('Invalid value for score. Allowed string values are "ba_ad", "rmse_ad".')

    def __make_inverse_matrix(self, X):
        X = column_stack(((ones(X.shape[0])), X))
        influence_matrix = X.T.dot(X) + eye(X.shape[1]).dot(1e-8)
        return linalg.inv(influence_matrix)

    def __find_leverages(self, X, inverse_influence_matrix):
        X = column_stack(((ones(X.shape[0])), X))
        return array([X[i, :].dot(inverse_influence_matrix).dot(X[i, :]) for i in range(X.shape[0])])

    def fit(self, X, y=None):
        """Learning is to find the inverse matrix for X and calculate the threshold.

        Parameters
        ----------
        X : array-like or sparse matrix, shape (n_samples, n_features)
            The input samples. Use ``dtype=np.float32`` for maximum
            efficiency.
        y : array-like, shape = [n_samples] or [n_samples, n_outputs]
            The target values (real numbers in regression).

        Returns
        -------
        self : object
        """
        # Check that X have correct shape
        X = check_array(X)
        self.inverse_influence_matrix = self.__make_inverse_matrix(X)
        if self.threshold == 'auto':
            self.threshold_value = 3 * (1 + X.shape[1]) / X.shape[0]
        elif self.threshold == 'cv':
            if y is None:
                raise ValueError("Y must be specified to find the optimal threshold.")
            y = check_array(y, accept_sparse='csc', ensure_2d=False, dtype=None)
            self.threshold_value = 0
            score_value = 0
            Y_pred, Y_true, AD = [], [], []
            cv = KFold(n_splits=5, random_state=1, shuffle=True)
            for train_index, test_index in cv.split(X):
                x_train = safe_indexing(X, train_index)
                x_test = safe_indexing(X, test_index)
                y_train = safe_indexing(y, train_index)
                y_test = safe_indexing(y, test_index)
                if self.reg_model is None:
                    reg_model = RandomForestRegressor(n_estimators=500, random_state=1).fit(x_train, y_train)
                else:
                    reg_model = clone(self.reg_model).fit(x_train, y_train)
                Y_pred.append(reg_model.predict(x_test))
                Y_true.append(y_test)
                ad_model = self.__make_inverse_matrix(x_train)
                AD.append(self.__find_leverages(x_test, ad_model))
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

    def predict_proba(self, X):
        """Predict the distances for X to center of the training set.

        Parameters
        ----------
        X : array-like or sparse matrix, shape (n_samples, n_features)
            The input samples. Internally, it will be converted to
            ``dtype=np.float32`` and if a sparse matrix is provided
            to a sparse ``csr_matrix``.

        Returns
        -------
        leverages: array of shape = [n_samples]
                   The objects distances to center of the training set.
        """
        # Check is fit had been called
        check_is_fitted(self, ['inverse_influence_matrix'])
        # Check that X have correct shape
        X = check_array(X)
        return self.__find_leverages(X, self.inverse_influence_matrix)

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
        check_is_fitted(self, ['inverse_influence_matrix', 'threshold_value'])
        # Check that X have correct shape
        X = check_array(X)
        return self.__find_leverages(X, self.inverse_influence_matrix) <= self.threshold_value


__all__ = ['Leverage']
