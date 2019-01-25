# -*- coding: utf-8 -*-
#
#  Copyright 2018 Assima Rakhimbekova <asima.astana@outlook.com>
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


class ReactionTypeControl():
    """
    Reaction Type Control (RTC) is performed using reaction signature.
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

    def __readFile(self, X):
        if self.env == 'all':
            return [str(~r) for r in X]
        else:
            data = []
            for r in X:
                cgr = ~r  # Condence Graph of Reaction
                cgr.reset_query_marks()  # reset hyb and neighbors marks to atoms
                center_atoms = cgr.get_center_atoms()  # Numbers atoms in reaction center in list
                aug_center = cgr.augmented_substructure(center_atoms, dante=False, deep=self.env,
                                                        as_view=True)  # get ubgraph with atoms and their neighbors
                data.append(format(aug_center, 'h'))  # String for graph reaction center
        return data

    def fit(self, X):
        """Fit structure-based AD.
        The training model  memorizes the unique set of reaction signature.
        Parameters
        ----------
        X : after read rdf file
        Returns
        -------
        self : object
            Returns self.
        """
        self.dataHash = set(self.__readFile(X))
        return self

    def predict(self, X):
        """
        Reaction is considered belonging to model’s AD if its reaction signature coincides with ones used in training set.
        :param X: after read rdf file.
        :return: list contains True (reaction in AD) and False (reaction residing outside AD).
        """
        state = [newHash in self.dataHash for newHash in self.__readFile(X)]
        return state


__all__ = ['ReactionTypeControl']
