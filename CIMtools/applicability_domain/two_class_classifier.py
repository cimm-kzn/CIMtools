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
from sklearn.base import BaseEstimator, clone
from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier
from sklearn.metrics import mean_squared_error
from sklearn.model_selection import cross_val_predict, KFold, GridSearchCV
from sklearn.utils import safe_indexing
from sklearn.utils.validation import check_array, check_is_fitted
from ..metrics.applicability_domain_metrics import balanced_accuracy_score_with_ad, rmse_score_with_ad


class TwoClassClassifiers(BaseEstimator):
    """
    Model learns to distinguish inliers from outliers.
    Objects with high prediction error in cross-validation (more than 3xRMSE) are considered outliers,
    while the rest are inliers.
    Two-class classification methods is trained to distinguish them, and provides the value of confidence that
    object belongs to inliers. The latter is used as a measure that object is in AD.
    In this case, Random Forest Classifier implemented in scikit-learn library is used.
    The method requires fitting of two hyperparameters: max_features and probability threshold P*.
    If the object’s predicted probability of belonging to the inliers is greater than P*,
    its prediction is considered reliable (within AD).
    Other hyperparameters of Random Forest Classifier were set to defaults,
    except number of decision trees in RF was set to 500.
    """
    def __init__(self, threshold='cv', score='ba_ad', reg_model=None, clf_model=None):
        self.threshold = threshold
        self.score = score
        self.reg_model = reg_model
        self.clf_model = clf_model
        if threshold != 'cv' and not isinstance(threshold, float):
            raise ValueError('Invalid value for threshold. Allowed string value is "cv".')
        if score not in ('ba_ad', 'rmse_ad'):
            raise ValueError('Invalid value for score. Allowed string values are "ba_ad", "rmse_ad".')

    def fit(self, X, y=None):
        """
        Model building and threshold searching
        During training, a model is built and a probability threshold is found by which the object is considered to
        belong to the applicability domain of the model.
        For this reason, in fit method we pass the following parameters:
        reg_model and clf_model. Reg_model is regression model, clf_model is classification model.

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
        if y is None:
            raise ValueError("Y must be specified to find the optimal threshold.")
        y = check_array(y, accept_sparse='csc', ensure_2d=False, dtype=None)
        # Check that X have correct shape
        X = check_array(X)

        if self.reg_model is None:
            reg_model = RandomForestRegressor(n_estimators=500, random_state=1)
        else:
            reg_model = clone(self.reg_model)
        if self.clf_model is None:
            clf_model = RandomForestClassifier(n_estimators=500, random_state=1)
        else:
            clf_model = clone(self.clf_model)

        cv = KFold(n_splits=5, random_state=1, shuffle=True)
        y_pred = cross_val_predict(reg_model, X, y, cv=cv)
        y_clf = abs(y_pred - y) <= 3 * sqrt(mean_squared_error(y, y_pred))
        self.AD_clf = clf_model.fit(X, y_clf)
        if self.threshold == 'cv':
            reg_model_int = clone(reg_model)
            clf_model_int = clone(clf_model)

            self.threshold_value = 0
            score_value = 0
            Y_pred, Y_true, AD = [], [], []
            for train_index, test_index in cv.split(X):
                x_train = safe_indexing(X, train_index)
                x_test = safe_indexing(X, test_index)
                y_train = safe_indexing(y, train_index)
                y_test = safe_indexing(y, test_index)
                y_train_clf = safe_indexing(y_clf, train_index)
                Y_pred.append(reg_model_int.fit(x_train, y_train).predict(x_test))
                Y_true.append(y_test)
                AD.append(clf_model_int.fit(x_train, y_train_clf).predict_proba(x_test)[:, 0])
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
        check_is_fitted(self, ['AD_clf', 'threshold_value'])
        # Check that X have correct shape
        X = check_array(X)
        return self.AD_clf.predict_proba(X)[:, 0] <= self.threshold_value


__all__ = ['TwoClassClassifiers']
