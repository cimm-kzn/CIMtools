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
from sklearn.base import BaseEstimator
from sklearn.utils.validation import check_array, check_is_fitted


class Box(BaseEstimator):
    """ This approach defines AD as a bounding block, which is an N-dimensional hypercube
    defined on the basis of the maximum and minimum values of each descriptor used to construct the model.
    If test compound is outside of hypercube it is outside of AD model.
    The method doesn’t have internal parameters, threshold.

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
        y : Ignored
            not used, present for API consistency by convention.

        Returns
        -------
        self : object
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


__all__ = ['Box']
