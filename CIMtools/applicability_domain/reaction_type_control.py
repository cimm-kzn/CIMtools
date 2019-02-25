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
from CGRtools.containers import ReactionContainer
from numpy import array
from sklearn.utils.validation import check_is_fitted
from ..utils import iter2array


class ReactionTypeControl:
    """Reaction Type Control (RTC) is performed using reaction signature.

    The signature includes both the reaction centre itself and its 1, 2, and so on the environment
    Since the reaction signature is not a very clear term, we considered the environment parameter as a hyper-parameter.
    Therefore, the method has one internal parameter. If the environment is 0,
    then the reaction signature is considered only the atoms at which the change occurs.
    If environment = 1, the first environment was included in the reaction signature,
    if environment = 2 - the second environment,
    and so on up to the whole reaction (env='all').
    In addition, by default, all atoms put a label on their hybridization.
    Reaction is considered belonging to model’s AD if its reaction signature coincides with ones used in training set.
    """
    def __init__(self, env=0):
        self.env = env

    def __get_signature(self, structure):
        if self.env == 'all':
            return str(~structure)
        else:
            cgr = ~structure  # Condence Graph of Reaction
            cgr.reset_query_marks()  # reset hyb and neighbors marks to atoms
            # get subgraph with atoms and their neighbors
            aug_center = cgr.augmented_substructure(cgr.center_atoms, deep=self.env)
            return format(aug_center, 'h')  # String for graph reaction center

    def fit(self, X):
        """Fit structure-based AD. The training model  memorizes the unique set of reaction signature.

        Parameters
        ----------
        X : after read rdf file

        Returns
        -------
        self : object
        """
        X = iter2array(X, dtype=ReactionContainer)
        self._train_signatures = {self.__get_signature(x) for x in X}
        return self

    def predict(self, X):
        """Reaction is considered belonging to model’s AD
        if its reaction signature coincides with ones used in training set.

        Parameters
        ----------
        X : after read rdf file

        Returns
        -------
        self : array contains True (reaction in AD) and False (reaction residing outside AD).
        """
        check_is_fitted(self, ['_train_signatures'])
        X = iter2array(X, dtype=ReactionContainer)
        return array([self.__get_signature(x) in self._train_signatures for x in X])


__all__ = ['ReactionTypeControl']
