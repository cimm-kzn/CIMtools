# -*- coding: utf-8 -*-
#
#
#  Copyright 2019 Assima Rakhimbekova <asima.astana@outlook.com>
#  This file is part of CIMtools.
#
#  CIMtools is free software; you can redistribute it and/or modify
#  it under the terms of the GNU Affero General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Affero General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, see <https://www.gnu.org/licenses/>.
#
from numpy import float64, hstack, sqrt, unique, zeros
from sklearn.base import clone
from sklearn.ensemble import *
from sklearn.ensemble.base import _partition_estimators
from sklearn.ensemble.forest import _accumulate_prediction
from sklearn.externals.joblib import *
from sklearn.model_selection import KFold
from sklearn.tree.tree import DecisionTreeRegressor
import threading
from sklearn.utils import safe_indexing
from sklearn.utils.validation import check_array, check_is_fitted
from ..metrics.applicability_domain_metrics import balanced_accuracy_score_with_ad, rmse_score_with_ad


class VariancesFromRandomForest(RandomForestRegressor):

    def __init__(self, score='ba_ad', threshold='cv', reg_model=None, n_estimators=100,
                 criterion="mse",
                 max_depth=None,
                 min_samples_split=2,
                 min_samples_leaf=1,
                 min_weight_fraction_leaf=0.,
                 max_features="auto",
                 max_leaf_nodes=None,
                 min_impurity_decrease=0.,
                 min_impurity_split=None,
                 bootstrap=True,
                 oob_score=False,
                 n_jobs=None,
                 random_state=None,
                 verbose=0,
                 warm_start=False):
        super(RandomForestRegressor, self).__init__(
            base_estimator=DecisionTreeRegressor(),
            n_estimators=n_estimators,
            estimator_params=("criterion", "max_depth", "min_samples_split",
                              "min_samples_leaf", "min_weight_fraction_leaf",
                              "max_features", "max_leaf_nodes",
                              "min_impurity_decrease", "min_impurity_split",
                              "random_state"),
            bootstrap=bootstrap,
            oob_score=oob_score,
            n_jobs=n_jobs,
            random_state=random_state,
            verbose=verbose,
            warm_start=warm_start)
        self.criterion = criterion
        self.max_depth = max_depth
        self.min_samples_split = min_samples_split
        self.min_samples_leaf = min_samples_leaf
        self.min_weight_fraction_leaf = min_weight_fraction_leaf
        self.max_features = max_features
        self.max_leaf_nodes = max_leaf_nodes
        self.min_impurity_decrease = min_impurity_decrease
        self.min_impurity_split = min_impurity_split
        self.criterion = criterion
        self.max_depth = max_depth
        self.min_samples_split = min_samples_split
        self.min_samples_leaf = min_samples_leaf
        self.min_weight_fraction_leaf = min_weight_fraction_leaf
        self.max_features = max_features
        self.max_leaf_nodes = max_leaf_nodes
        self.min_impurity_decrease = min_impurity_decrease
        self.min_impurity_split = min_impurity_split
        self.score = score
        self.threshold = threshold
        self.reg_model = reg_model
        if threshold is not'cv' and not isinstance(threshold, float):
            raise ValueError('Invalid value for threshold. Allowed string value is "cv".')
        if score not in ('ba_ad', 'rmse_ad'):
            raise ValueError('Invalid value for score. Allowed string values are "ba_ad", "rmse_ad".')

    def _accumulate_prediction2(self, predict, X, out, lock):
        prediction = predict(X, check_input=False)
        with lock:
            if len(out) == 1:
                out[0] += (prediction ** 2)
            else:
                for i in range(len(out)):
                    out[i] += (prediction[i] ** 2)

    def _variance_of_values(self, X, estimators_):
        """Predict regression target for X.

        The predicted regression target of an input sample is computed as the
        mean predicted regression targets of the trees in the forest.

        Parameters
        ----------
        X : array-like or sparse matrix of shape = [n_samples, n_features]
            The input samples. Internally, its dtype will be converted to
            ``dtype=np.float32``. If a sparse matrix is provided, it will be
            converted into a sparse ``csr_matrix``.

        Returns
        -------
        variances : array of shape = [n_samples] or [n_samples, n_outputs]
                    The predicted variance of values.
        """
        check_is_fitted(self, 'estimators_')
        # Check data
        X = self._validate_X_predict(X)

        # Assign chunk of trees to jobs
        n_jobs, _, _ = _partition_estimators(self.n_estimators, self.n_jobs)

        # avoid storing the output of every estimator by summing them here
        if self.n_outputs_ > 1:
            y_hat = zeros((X.shape[0], self.n_outputs_), dtype=float64)
            y_hat_new = zeros((X.shape[0], self.n_outputs_), dtype=float64)
        else:
            y_hat = zeros((X.shape[0]), dtype=float64)
            y_hat_new = zeros((X.shape[0]), dtype=float64)

        # Parallel loop
        lock = threading.Lock()
        Parallel(n_jobs=n_jobs, verbose=self.verbose, backend="threading")(
            delayed(_accumulate_prediction)(e.predict, X, [y_hat], lock)
            for e in estimators_)

        # Parallel loop 2
        lock2 = threading.Lock()
        Parallel(n_jobs=n_jobs, verbose=self.verbose, backend="threading")(
            delayed(self._accumulate_prediction2)(e.predict, X, [y_hat_new], lock2)
            for e in estimators_)

        y_hat /= len(estimators_)
        y_hat_new /= len(estimators_)
        return sqrt(y_hat_new - (y_hat ** 2))

    def fit(self, X, y=None):
        self.estimators_ = super(RandomForestRegressor, self).fit(X, y)
        if self.threshold == 'cv':
            if y is None:
                raise ValueError("Y must be specified to find the optimal threshold.")
            y = check_array(y, accept_sparse='csc', ensure_2d=False, dtype=None)
            self.threshold_value = 0
            score = 0
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
                ad_model = super(RandomForestRegressor, self).fit(x_train, y_train)
                AD.append(self._variance_of_values(x_test, ad_model))
            AD_ = unique(hstack(AD))
            for z in AD_:
                AD_new = hstack(AD) <= z
                if self.score == 'ba_ad':
                    val = balanced_accuracy_score_with_ad(Y_true=hstack(Y_true), Y_pred=hstack(Y_pred), AD=AD_new)
                elif self.score == 'rmse_ad':
                    val = rmse_score_with_ad(Y_true=hstack(Y_true), Y_pred=hstack(Y_pred), AD=AD_new)
                if val >= score:
                    score = val
                    self.threshold_value = z
        else:
            self.threshold_value = self.threshold
        return self

    def predict(self, X):
        check_is_fitted(self, 'estimators_')
        # Check data
        X = self._validate_X_predict(X)
        y = self._variance_of_values(X, self.estimators_)
        return y <= self.threshold_value


__all__ = ['VariancesFromRandomForest']
