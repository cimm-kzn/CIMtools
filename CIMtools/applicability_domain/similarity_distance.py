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
from numpy import hstack, mean, sqrt, var, unique
from sklearn.base import BaseEstimator, clone
from sklearn.ensemble import RandomForestRegressor
from sklearn.neighbors import BallTree
from sklearn.model_selection import KFold
from sklearn.utils import safe_indexing
from sklearn.utils.validation import check_array, check_is_fitted
from ..metrics.applicability_domain_metrics import balanced_accuracy_score_with_ad, rmse_score_with_ad


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

    NB! To the nearest first neighbor

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

    threshold : string or float
        It needs to compare the distance values with threshold. If the calculated distance values of test set compounds
        are not within the threshold set by the training set molecules, then the prediction of these
        compounds are considered to be unreliable.

        - If auto, threshold calculated like Dc = Zσ + <y>, where <y> is the average and σ is the standard deviation of
            the Euclidean distances of the k nearest neighbors of each compound in the training set and
            Z is an empirical parameter to control the significance level, with the default value of 0.5.
        - If 'cv', threshold in the approach is optimized in course internal cross-validation procedure
            by maximize our metric.
        - IF float, threshold will be this value

    score : string
        A metric is required to find a threshold.

        - If score is 'ba_ad' is calculated balanced accuracy. The true inliers and outliers are those for
            which the difference in the prediction error is less than 3 RMSE
        - If score is 'rmse_ba' is calculated Root Mean Squared Error of model with AD. Sahigata and etc proposed [1]
            to use difference between root mean squared error outliers and inliers (RMSE_AD), which shows
            what is predicted better: objects outside AD or objects inside and outside AD.
            The metric characterizes how accurate the model becomes. By inliers, we mean objects inside AD,
            and by outliers, objects outside AD.

    reg_model : None or estimator
        It needs for finding threshold

    ----
    [1] Sahigara F., Mansouri K., Ballabio D., Mauri A., Consonni V. Todeschini R. Comparison of Different Approaches
    to Define the Applicability Domain of QSAR Models.  Molecules, 2012, vol. 17, pp. 4791-4810.
    doi: 10.3390/molecules17054791.
    """
    def __init__(self, leaf_size=40, metric='minkowski', score='ba_ad', threshold='auto', reg_model=None):
        self.leaf_size = leaf_size
        self.metric = metric
        self.score = score
        self.threshold = threshold
        self.reg_model = reg_model
        if threshold not in ('auto','cv') and not isinstance(threshold, float):
            raise ValueError('Invalid value for threshold. Allowed string values are "auto", "cv".')
        if score not in ('ba_ad', 'rmse_ad'):
            raise ValueError('Invalid value for score. Allowed string values are "ba_ad", "rmse_ad".')

    def fit(self, X, y=None):
        """Fit distance-based AD.

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
        self.tree = BallTree(X, leaf_size=self.leaf_size, metric=self.metric)
        dist_train = self.tree.query(X, k=2)[0]
        if self.threshold == 'auto':
            self.threshold_value = 0.5 * sqrt(var(dist_train[:, 1])) + mean(dist_train[:, 1])
        elif self.threshold == 'cv':
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
                data_test = safe_indexing(dist_train[:, 1], test_index)
                if self.reg_model is None:
                    reg_model = RandomForestRegressor(n_estimators=500, random_state=1).fit(x_train, y_train)
                else:
                    reg_model = clone(self.reg_model).fit(x_train, y_train)
                Y_pred.append(reg_model.predict(x_test))
                Y_true.append(y_test)
                AD.append(data_test)
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

    def predict_proba(self, X):
        """Returns the value of the nearest neighbor from the training set.

        Parameters
        ----------
         X : array-like or sparse matrix, shape (n_samples, n_features)
            The input samples. Internally, it will be converted to
            ``dtype=np.float32`` and if a sparse matrix is provided
            to a sparse ``csr_matrix``.

        Returns
        -------
        y : array, shape (n_samples,)
        """
        # Check is fit had been called
        check_is_fitted(self, ['tree'])
        # Check data
        X = check_array(X)
        return self.tree.query(X)[0].flatten()

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
        y : array, shape (n_samples,)
            For each observations, tells whether or not (True or False) it should
            be considered as an inlier according to the fitted model.
        """
        # Check is fit had been called
        check_is_fitted(self, ['tree'])
        # Check data
        X = check_array(X)
        return self.tree.query(X)[0].flatten() <= self.threshold_value


__all__ = ['SimilarityDistance']
