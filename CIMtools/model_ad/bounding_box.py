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
from sklearn.base import BaseEstimator
from sklearn.utils.validation import check_array, check_is_fitted


class Box(BaseEstimator):
    """ This approach defines AD as a bounding block, which is an N-dimensional hypercube
    defined on the basis of the maximum and minimum values of each descriptor used to construct the model.
    If test compound is outside of hypercube it is outside of AD model.
    The method doesn’t have internal parameters, threshold.
    Пример
    ------
    reactions = RDFread(file).read()
    for train_index_ext, test_index_ext in kf.split(reactions):  # external set
        reactions_train = safe_indexing(reactions, train_index_ext)
        reactions_test = safe_indexing(reactions, test_index_ext)
        data = initial_data(reactions_train, constant_name, reactions_test, constant_name)
        X_train = data[0][0] # np.array
        Y_train = data[0][1] # list
        y_train_gpr = data[0][2][:, 0] # np.array
        X_test = data[1][0]  # np.array
        Y_test = data[1][1]  # list
        y_test_gpr = data[1][2][:, 0]  # np.array
        ... строим модель регрисионную ...

        # **************************************** REACTION TYPE CONTROL DEEP = 1 **************************************
        RTC_1 = ReactionTypeControl(env=1, hyb=True).fit(reactions_train)
        AD_RTC_deep1 = RTC_1.predict(reactions_test)
        res(AD_res, 0, AD_RTC_deep1)
        # ******************************************** FRAGMENT CONTROL ************************************************
        FC = Pipeline([('cgr', CGR()), ('frg', Fragmentor(version='2017.x', max_length=4, useformalcharge=True,
                                                          return_domain=True))]).fit(reactions_train)
        FC_AD = FC.transform(reactions_test)
        AD_FC_results = FC_AD.iloc[:, FC_AD.shape[1] - 1]
        AD_FC = AD_FC_results.get_values()
        res(AD_res, 1, AD_FC) #
        #********************************************** BOUNDING BOX ***************************************************
        BB = Box().fit(X_train)
        AD_BB = BB.predict(X_test)
        res(AD_res, 6, AD_BB, FC=AD_FC, RTC_deep_1=AD_RTC_deep1)

    Parameters
    ----------
    Bounding Box doesn't have parameters.
    """

    def __init__(self):
        pass

    def fit(self, X, y=None):
        """Find min and max values of every feature.
        Parameters
        ----------
        X : {array-like, sparse matrix}, shape (n_samples, n_features)
            The training input samples.
        Returns
        -------
        self : object
            Returns self.
        """
        # Check that X have correct shape
        X = check_array(X)

        self._x_min = X.min(axis=0) # axis=0 will find the minimum values ​​by columns (for each feature)
        self._x_max = X.max(axis=0) # axis=0 will find the minimum values ​​by columns (for each feature)
        return self

    def predict(self, X):
        """ Predict if a particular sample is an outlier or not.

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
        # Check is fit had been called
        check_is_fitted(self, ['_x_min', '_x_max'])

        # Input validation
        X = check_array(X)
        return ((X - self._x_min).min(axis=1) >= 0) & ((self._x_max - X).min(axis=1) >= 0)
