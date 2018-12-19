# -*- coding: utf-8 -*-
#
#  Copyright 2018 Tagir Akhmetshin <tagirshin@gmail.com>
#  Copyright 2018 Ramil Nugmanov <stsouko@live.ru>
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
from collections import defaultdict
from numpy import array
from sklearn.model_selection import BaseCrossValidator
from sklearn.utils.validation import indexable


class LeaveOneGroupOut(BaseCrossValidator):
    """Leave-One-Group-Out cross-validator

    Provides train/test indices to split data in train/test sets. Each
    reactions with the same condition is used once as a test set (singleton)
    while the remaining reactions form the training set. Test set includes
    only reactions with transformations that appeared in other reactions.
    """
    def get_n_splits(self, X=None, y=None, groups=None):
        """Returns the number of splitting iterations in the cross-validator
        Parameters
        ----------
        X : object
            Always ignored, exists for compatibility.
            ``np.zeros(n_samples)`` may be used as a placeholder.
        y : object
            Always ignored, exists for compatibility.
            ``np.zeros(n_samples)`` may be used as a placeholder.
        groups : array-like, with shape (n_samples,)
            Group labels for the samples used while splitting the dataset into
            train/test set.
        Returns
        -------
        n_splits : int
            Returns the number of splitting iterations in the cross-validator.
        """
        return len(set(groups))

    def split(self, X, y=None, groups=None):
        """Generate indices to split data into training and test set.
        Parameters
        ----------
        X : array-like, of length n_samples
            Training data, includes reaction's containers
        y : array-like, of length n_samples
            The target variable for supervised learning problems.
        groups : array-like, with shape (n_samples,)
            Group labels for the samples used while splitting the dataset into
            train/test set.
        Yields
        ------
        train : ndarray
            The training set indices for that split.
        test : ndarray
            The testing set indices for that split.
        """
        X, y, groups = indexable(X, y, groups)
        cgrs = [~r for r in X]

        structure_condition = defaultdict(set)

        for structure, condition in zip(cgrs, groups):
            structure_condition[structure].add(condition)

        train_data = defaultdict(list)
        test_data = []

        for n, (structure, condition) in enumerate(zip(cgrs, groups)):
            train_data[condition].append(n)
            if len(structure_condition[structure]) > 1:
                test_data.append(n)

        for condition, indexes in train_data.items():
            test_index = [index for index in indexes if index in test_data]
            if test_index:
                train_index = [i for cond, ind in train_data.items() if cond != condition for i in ind]
                yield array(train_index), array(test_index)


__all__ = ['LeaveOneGroupOut']
