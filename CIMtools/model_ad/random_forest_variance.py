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
import numpy as np
from sklearn.ensemble import *
from sklearn.ensemble.forest import _accumulate_prediction
from sklearn.ensemble.forest import check_is_fitted
from sklearn.externals.joblib import *
from sklearn.ensemble.base import _partition_estimators
import threading
from sklearn.tree.tree import DecisionTreeRegressor


def _accumulate_prediction2(predict, X, out, lock):
    prediction = predict(X, check_input=False)
    with lock:
        if len(out) == 1:
            out[0] += (prediction ** 2)
        else:
            for i in range(len(out)):
                out[i] += (prediction[i] ** 2)


class RandomForestRegressor2(RandomForestRegressor):
    """
    Модели AD, основанные на согласие множества деревьев. для предсказания свойств химических реакций использовали
    случайный лес (500 деревьев), согласие каждого дерева, то есть дисперсия предсказанных значений каждого дерева
    также может рассматриваться как способ оценки AD. Гипер-параметром такой модели является также пороговое значение,
    если для объекта дисперсия предсказанных значений меньше порогового, то объект принадлежит AD. RFR_var дисперсия
    предсказанных значений полученных из 500 деревьев.

    # *********************************** RFR_VAR ******************************************************************
        est_var = RandomForestRegressor2(random_state=seed, n_estimators=500, max_features=est.best_params_['max_features'],
                                         n_jobs=4).fit(X_train, Y_train)
        AD_est_var_values = est_var.predict_proba(X_test)
        min_h_param_RFR_VAR = threshold(ad='rfr_var', X=X_train, y=Y_pr_ts, metric=score, reg_model=est.best_estimator_,
                                        ad_model=est_var) # для нахождения отсечки
        AD_var_RF_itog = AD_est_var_values <= min_h_param_RFR_VAR['z']

    """
    def __init__(self, n_estimators=100,
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

    def predict_proba(self, X):
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
        y : array of shape = [n_samples] or [n_samples, n_outputs]
            The predicted values.
        """
        check_is_fitted(self, 'estimators_')
        # Check data
        X = self._validate_X_predict(X)

        # Assign chunk of trees to jobs
        n_jobs, _, _ = _partition_estimators(self.n_estimators, self.n_jobs)

        # avoid storing the output of every estimator by summing them here
        if self.n_outputs_ > 1:
            y_hat = np.zeros((X.shape[0], self.n_outputs_), dtype=np.float64)
            y_hat_new = np.zeros((X.shape[0], self.n_outputs_), dtype=np.float64)
        else:
            y_hat = np.zeros((X.shape[0]), dtype=np.float64)
            y_hat_new = np.zeros((X.shape[0]), dtype=np.float64)

        # Parallel loop
        lock = threading.Lock()
        Parallel(n_jobs=n_jobs, verbose=self.verbose, backend="threading")(
            delayed(_accumulate_prediction)(e.predict, X, [y_hat], lock)
            for e in self.estimators_)

        # Parallel loop 2
        lock2 = threading.Lock()
        Parallel(n_jobs=n_jobs, verbose=self.verbose, backend="threading")(
            delayed(_accumulate_prediction2)(e.predict, X, [y_hat_new], lock2)
            for e in self.estimators_)

        y_hat /= len(self.estimators_)
        y_hat_new /= len(self.estimators_)
        var = np.sqrt(y_hat_new - (y_hat ** 2))
        return var
