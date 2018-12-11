# -*- coding: utf-8 -*-
#
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
from CGRtools import CGRpreparer
from CGRtools.containers import ReactionContainer
from sklearn.base import BaseEstimator
from ..common import iter2array, nested_iter_to_2d_array, TransformerMixin


class AtomMarkerCGR(BaseEstimator, TransformerMixin):
    def __init__(self, templates, only_first=True):
        """
        :param templates: list of tuples of Query and marks dict
        :param only_first: return only first match
        """
        self.templates = tuple(templates)
        self.only_first = only_first
        self.__cgr = CGRpreparer()

    def __getstate__(self):
        return {k: v for k, v in super().__getstate__().items() if not k.startswith('_AtomMarkerCGR__')}

    def __setstate__(self, state):
        super().__setstate__(state)
        self.__cgr = CGRpreparer()

    def transform(self, x):
        x = super().transform(x)

        if self.only_first:
            return iter2array((self.__prepare(r) or None for r in x), allow_none=True)
        return nested_iter_to_2d_array((self.__prepare(r) for r in x), allow_none=True)

    def __prepare(self, structure):
        cgr = self.__cgr.compose(structure)
        result = []
        for query, marks in self.templates:
            for mapping in query.get_substructure_mapping(cgr, -1):
                s = structure.copy()
                for k, v in marks.items():
                    atom = mapping[k]
                    r = next((x for x in s.reagents if atom in x), None)
                    p = next((x for x in s.products if atom in x), None)
                    if r:
                        r.atom(atom).mark = v
                    if p:
                        p.atom(atom).mark = v
                if self.only_first:
                    return s
                result.append(s)

        return result

    _dtype = ReactionContainer


__all__ = ['AtomMarkerCGR']
