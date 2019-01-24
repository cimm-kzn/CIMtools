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
from sklearn.neighbors import BallTree
from sklearn.utils.validation import check_array, check_is_fitted
from sklearn.base import BaseEstimator
from sklearn.ensemble import RandomForestRegressor
from CIMtools.applicability_domain.domain_selection.threshold_functions import threshold


class SimilarityDistance(BaseEstimator):
    """ Distance-based method for  defining applicability domain (AD).
    In the case of non-linear kNN QSPR method, since the models are based on chemical similarity calculations,
    a large similarity distance could signal query compounds too dissimilar to the training set compounds.
    This approach is based on providing similarity measure for a new chemical with respect to the compounds within
    the training space. The similarity is identified by finding the distance of a query chemical from the nearest
    training compound or its distances from k nearest neighbors in the training set.
    If the calculated distance values of test set compounds are not within the user-defined threshold set by
    the training set molecules, then the prediction of these compounds are considered to be unreliable.
    Commonly threshold calculated like Dc=Zσ + <y>, where <y> is the average and σ is the standard deviation of the
    Euclidean distances of the k nearest neighbors of each compound in the training set and Z is an empirical parameter
    to control the significance level, with the default value of 0.5.

    Drawback of method is lack of strict rules in literature towards defining the thresholds can lead to ambiguous
    results. We propose a variation of finding threshold. Threshold in the approach was optimized in course internal
    cross-validation procedure by maximize our metric.

    NB! До ближайщего первого соседа

    Read more in the article.

    Parameters
    ----------
    leaf_size : positive integer (default = 40)
        Number of points at which to switch to brute-force. Changing leaf_size will not affect the results of a query,
        but can significantly impact the speed of a query and the memory required to store the constructed tree.
        The amount of memory needed to store the tree scales as approximately n_samples / leaf_size. For a specified
        leaf_size, a leaf node is guaranteed to satisfy leaf_size <= n_points <= 2 * leaf_size, except in the case that
        n_samples < leaf_size.

    metric : string or DistanceMetric object
        The distance metric to use for the tree. Default=’minkowski’ with p=2 (that is, a euclidean metric).
        See the documentation of the DistanceMetric class for a list of available metrics. ball_tree.valid_metrics gives
         a list of the metrics which are valid for BallTree.

    threshold
        It needs to compare the distance values with threshold. If the calculated distance values of test set compounds
        are not within the user-defined threshold set by the training set molecules, then the prediction of these
        compounds are considered to be unreliable.

        - If None, threshold calculated like Dc=Zσ + <y>, where <y> is the average and σ is the standard deviation of
            the Euclidean distances of the k nearest neighbors of each compound in the training set and Z is an empirical
            parameter to control the significance level, with the default value of 0.5.
        - If 'threshold_Cv', threshold in the approach was optimized in course internal cross-validation procedure
            by maximize our metric.
    --------
    В итоге сделала так, чтобы два варианта сразу рассчитывались

    Пример
    ------
        # #******************************************** Similatiry Distance ********************************************
    SIM = SimilarityDistance(ad='simdist', mscore=score, reg_model=est.best_estimator_).fit(X_train, Y_pr_ts)
    SIM_pred = SIM.predict(X_test)
    res(AD_res, 18, SIM_pred[:, 0], FC=AD_FC, RTC_deep_1=AD_RTC_deep1)
    res(AD_res, 21, SIM_pred[:, 1], FC=AD_FC, RTC_deep_1=AD_RTC_deep1)
    """
    def __init__(self, leaf_size=40, metric='minkowski',  ad='simdist', mscore='ba',
                 reg_model=RandomForestRegressor(n_estimators=500, random_state=1)):
        self.leaf_size = leaf_size
        self.mscore = mscore # у меня везде написано metric, но тут еще один параметр также называется. Упсссс.
        # По умолчанию 'ba'
        self.ad = ad
        self.metric = metric
        self.reg_model = reg_model #  необходима для нахождения отсечки, если пользователь не задает регресиооную модель,
        # то по умолчанию будет RandomForestRegressor(n_estimators=500, random_state=1). Можно также добавить внутри
        # gridSearchCV если надо


    def fit(self, X, y):
        """Fit distance-based AD.
        В качестве обучения - это расчет значений расстояний каждого объекта до соседа

        Parameters
        ----------
        X : array-like or sparse matrix, shape (n_samples, n_features)
            The input samples. Use ``dtype=np.float32`` for maximum
            efficiency.

        Returns
        -------
        self : object
            Returns self.
        """
        # Check data
        X = check_array(X)
        y = check_array(y, accept_sparse='csc', ensure_2d=False, dtype=None) # нужен не для обучения, а для нахождения отсечки
        self.tree = BallTree(X, leaf_size=self.leaf_size, metric=self.metric)
        self.dist_train, ind_train = self.tree.query(X, k=2) # k = 2 два соседа, 1 - сам с собой, 2 - ближайший первый сосед
        self.average = np.mean(self.dist_train[:, 1]) # dist_train[:, 1] = нули, так как расстояние до себя самого же 0
        self.standard_deviation = np.sqrt(np.var(self.dist_train[:, 1]))
        self.simdist_threshold = threshold(ad=self.ad, X=X, y=y, metric=self.mscore, reg_model=self.reg_model,
                                            ad_model=None, envs=None, y_clf=None, data=self.dist_train[:, 1])
        return self

    def predict_proba(self, X):
        """Predict if a particular sample is an outlier or not.

        Parameters
        ----------
         X : array-like or sparse matrix, shape (n_samples, n_features)
           The input samples. Internally, it will be converted to
           ``dtype=np.float32`` and if a sparse matrix is provided
           to a sparse ``csr_matrix``.

       Returns
        -------
       is_inlier : array, shape (n_samples,)
           For each observations, tells whether or not (True or False) it should
           be considered as an inlier according to the fitted model.
        """
        # Check data
        X = check_array(X)

        # Check is fit had been called
        check_is_fitted(self, ['tree', 'simdist_threshold'])

        dist_test, ind_test = self.tree.query(X, k=1) # вычисляется расстояние до ближай
        return dist_test[:, 0]

    def predict(self, X):
        """Predict if a particular sample is an outlier or not.

        Parameters
        ----------
         X : array-like or sparse matrix, shape (n_samples, n_features)
           The input samples. Internally, it will be converted to
           ``dtype=np.float32`` and if a sparse matrix is provided
           to a sparse ``csr_matrix``.

       Returns
        -------
       is_inlier : array, shape (n_samples,)
           For each observations, tells whether or not (True or False) it should
           be considered as an inlier according to the fitted model.
        """
        # Check data
        X = check_array(X)

        # Check is fit had been called
        check_is_fitted(self, ['tree', 'simdist_threshold'])

        dist_test, ind_test = self.tree.query(X, k=1) # вычисляется расстояние до ближай
        AD_SIM = dist_test[:, 0] <= self.simdist_threshold['z']
        AD_similarity_dist = dist_test[:, 0] <= (0.5 * self.standard_deviation + self.average)
        return np.column_stack((AD_SIM, AD_similarity_dist))