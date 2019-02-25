# -*- coding: utf-8 -*-
#
#  Copyright 2018, 2019 Ramil Nugmanov <stsouko@live.ru>
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
from itertools import tee, chain
from sklearn.base import TransformerMixin
from .utils import iter2array


class CIMtoolsTransformerMixin(TransformerMixin):
    def fit(self, x, y=None):
        """Do nothing and return the estimator unchanged

        This method is just there to implement the usual API and hence work in pipelines.
        """
        if self._dtype is not None:
            iter2array(x, dtype=self._dtype)
        else:
            iter2array(x)
        return self

    def transform(self, x):
        if self._dtype is not None:
            return iter2array(x, dtype=self._dtype)
        return iter2array(x)

    _dtype = None


def reaction_support(_class):
    class ReactionSupport(_class):
        def transform(self, x):
            assert all(isinstance(s, ReactionContainer) for s in x), 'invalid dtype, olny ReactionContainers acceptable'

            shifts = {}
            mols = []
            for i in ('reactants', 'products'):
                sh = shifts[i] = [len(mols)]
                for s in x:
                    si = s[i]
                    sh.append(len(si) + sh[-1])
                    mols.extend(si)

            transformed = super().transform(mols)
            assert len(transformed) == len(mols), 'unexpected transformed molecules amount'

            out = []
            for s, r, p in zip(x, (transformed[y: z] for y, z in self.__pairwise(shifts['reactants'])),
                                  (transformed[y: z] for y, z in self.__pairwise(shifts['products']))):
                if any(i is None for i in chain(r, p)):
                    out.append(None)
                else:
                    out.append(ReactionContainer(r, p, meta=s.meta))
            return iter2array(out, allow_none=True)

        @staticmethod
        def __pairwise(iterable):
            """s -> (s0,s1), (s1,s2), (s2, s3), ..."""
            a, b = tee(iterable)
            next(b, None)
            return zip(a, b)

    return ReactionSupport
