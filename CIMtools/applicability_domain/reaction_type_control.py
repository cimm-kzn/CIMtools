# -*- coding: utf-8 -*-
#
#  Copyright 2019 Assima Rakhimbekova <asima.astana@outlook.com>
#  Copyright 2019, 2020 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from sklearn.base import BaseEstimator, ClassifierMixin
from sklearn.utils.validation import check_is_fitted
from ..utils import iter2array


class ReactionTypeControl(BaseEstimator, ClassifierMixin):
    """Reaction Type Control (RTC) is performed using reaction signature.

    The signature includes both the reaction centre itself and its nearest environment up to {env}
    Since the reaction signature is not a very clear term, we considered the environment parameter as a hyper-parameter.
    Therefore, the method has one internal parameter. If the environment is 0,
    then the reaction signature considers only the atoms at which the change occurs.
    If environment = 1, the first circle neighbours included in the reaction signature,
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
            # get subgraph with atoms and their neighbors
            aug_center = cgr.augmented_substructure(cgr.center_atoms, deep=self.env, as_query=True)
            # remove neighbors marks
            sn = aug_center._neighbors
            pn = aug_center._p_neighbors
            for n in aug_center:
                sn[n] = pn[n] = ()
            return str(aug_center)  # String for graph reaction center

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
        a : array contains True (reaction in AD) and False (reaction residing outside AD).
        """
        check_is_fitted(self, ['_train_signatures'])
        X = iter2array(X, dtype=ReactionContainer)
        return array([self.__get_signature(x) in self._train_signatures for x in X])


__all__ = ['ReactionTypeControl']
