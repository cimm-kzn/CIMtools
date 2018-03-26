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
from CGRtools import CGRreactor
from CGRtools.containers import ReactionContainer
from sklearn.base import BaseEstimator, TransformerMixin
from ..common import iter2array


class StandardizeCGR(BaseEstimator, TransformerMixin):
    def __init__(self, templates, balance_groups=False, extralabels=False, isotope=False, element=True, stereo=False):
        """
        CGR standardization and reaction balancing

        :param templates: CGRTemplates. rules for graph modifications. possible be False
        :param balance_groups: if True: for unbalanced reactions contains multiple attached functional groups in
            products and one of them described in reagents - will be restored information about all equal groups.
            for example:

                R + B1-X-> B'1-R'-B'2 + X'

            where B' is transformed B, R and X same.
            we know what B'1 and B'2 is equal and B'1 is transformed B1 =>
            this groups most likely appeared from a single reagent. we can add copy of B-X to reagents.
            results will be:

                R + B1-X1 + B2-X2 -> B'1-R'-B'2 + X'1 + X'2

        :param extralabels: see CGRreactor init
        :param isotope: see CGRreactor init
        :param element: see CGRreactor init
        :param stereo: see CGRreactor init
        """
        self.templates = templates
        self.balance_groups = balance_groups
        self.extralabels = extralabels
        self.isotope = isotope
        self.element = element
        self.stereo = stereo

        assert templates or balance_groups, 'invalid params. need balance_groups or/and templates'

    def get_params(self, *args, **kwargs):
        params = super().get_params(*args, **kwargs)
        params['templates'] = [x.pickle() for x in params['templates']]
        return params

    def set_params(self, templates, **params):
        return super().set_params(templates=[ReactionContainer.unpickle(x) for x in templates], **params)

    def transform(self, x, y=None):
        if self.__reactor is None:
            self.__reactor = CGRreactor(extralabels=self.extralabels, isotope=self.isotope, element=self.element,
                                        stereo=self.stereo)
            if self.__searcher is None and self.templates:
                self.__searcher = self.__reactor.get_template_searcher(self.__reactor.prepare_templates(self.templates))

        return iter2array(self.__prepare(g) if g is not None else None for g in iter2array(x))

    def __prepare(self, g):
        if self.balance_groups:
            g = self.__reactor.clone_subgraphs(g)

        cond = self.__searcher is not None
        report = []
        while cond:
            searcher = self.__searcher(g)
            first_match = next(searcher, None)
            if not first_match:
                if report:
                    g.graph.setdefault('CGR_REPORT', []).extend(report)
                break

            g = self.__reactor.patcher(g, first_match.patch)
            if 'CGR_TEMPLATE' in first_match.meta:
                report.append(first_match.meta['CGR_TEMPLATE'])

            for match in searcher:
                g = self.__reactor.patcher(g, match.patch)
                if 'CGR_TEMPLATE' in match.meta:
                    report.append(match.meta['CGR_TEMPLATE'])
        return g

    __searcher = __reactor = None
