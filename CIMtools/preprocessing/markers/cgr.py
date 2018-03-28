# -*- coding: utf-8 -*-
#
#  Copyright 2018 Ramil Nugmanov <stsouko@live.ru>
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
#  You should have received a copy of the GNU Affero General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#
from CGRtools.containers import ReactionContainer
from CGRtools.preparer import CGRpreparer
from CGRtools.reactor import CGRreactor
from sklearn.base import BaseEstimator
from ..common import iter2array, nested_iter_to_2d_array, TransformerMixin


class AtomMarkerCGR(BaseEstimator, TransformerMixin):
    def __init__(self, templates, extralabels=False, isotope=False, element=True, stereo=False, only_first=True):
        self.templates = templates
        self.element = element
        self.isotope = isotope
        self.stereo = stereo
        self.extralabels = extralabels
        self.only_first = only_first
        self.__load()

    def __getstate__(self):
        return {k: v for k, v in super().__getstate__().items() if not k.startswith('_AtomMarkerCGR__')}

    def __setstate__(self, state):
        super().__setstate__(state)
        self.__load()

    def __load(self):
        self.__cgr = CGRpreparer(extralabels=self.extralabels)
        self.__react = CGRreactor(stereo=self.stereo, extralabels=self.extralabels, isotope=self.isotope,
                                  element=self.element)
        templates = self.__react.prepare_templates(self.templates)
        markers = len([x for _, x in templates[0].patch.nodes(data='mark') if x != '0'])
        assert markers, 'marks not found in templates'

        self.__markers = markers
        self.__templates = self.__react.get_template_searcher(templates)

    def transform(self, x):
        x = super().transform(x)

        if self.only_first:
            return iter2array((self.__prepare(r) or None for r in x), allow_none=True)
        return nested_iter_to_2d_array((self.__prepare(r) for r in x), allow_none=True)

    def get_count(self):
        return self.__markers

    def __prepare(self, structure):
        cgr = self.__cgr.condense(structure)
        result = []
        for match in self.__templates(cgr):
            s = structure.copy()
            for atom, a_mark in match.patch.nodes(data='mark'):
                r = next((x for x in s.reagents if atom in x), None)
                p = next((x for x in s.products if atom in x), None)
                if r:
                    r.node[atom]['mark'] = a_mark
                if p:
                    p.node[atom]['mark'] = a_mark
            if self.only_first:
                return s

            result.append(s)

        return result

    __cgr = __react = __markers = __templates = None
    _dtype = ReactionContainer
