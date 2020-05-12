# -*- coding: utf-8 -*-
#
#  Copyright 2018-2020 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from itertools import tee
from sklearn.base import BaseEstimator, TransformerMixin
from .utils import iter2array


class CIMtoolsTransformerMixin(TransformerMixin, BaseEstimator):
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


def reaction_support(_class, module=None):
    """
    Class Factory for transformers without reactions support.

    :param module: Current module name. required for pickle. By default _class module used.
    """
    class ReactionSupported(_class):
        def fit(self, x, y=None, **fit_params):
            return self.__run(False, x, y=y, **fit_params)

        def transform(self, x):
            return self.__run(True, x)

        def fit_transform(self, x, y=None, **fit_params):
            self.fit(x, y, **fit_params)
            return self.transform(x)

        def __run(self, transform, x, **kwargs):
            x = iter2array(x, dtype=ReactionContainer)

            mols = []
            r_shifts = [0]
            for s in x:
                ms = s.reactants
                r_shifts.append(len(ms) + len(mols))
                mols.extend(ms)

            p_shifts = [len(mols)]
            for s in x:
                ms = s.products
                p_shifts.append(len(ms) + len(mols))
                mols.extend(ms)

            g_shifts = [len(mols)]
            for s in x:
                ms = s.reagents
                g_shifts.append(len(ms) + len(mols))
                mols.extend(ms)

            if transform:
                transformed = super().transform(mols)
                if len(transformed) != len(mols):
                    raise ValueError('unexpected transformed molecules amount')

                return [(r, p, g) for r, p, g in
                        zip((transformed[y: z] for y, z in self.__pairwise(r_shifts)),
                            (transformed[y: z] for y, z in self.__pairwise(p_shifts)),
                            (transformed[y: z] for y, z in self.__pairwise(g_shifts)))]
            else:
                return super().fit(mols, **kwargs)

        @staticmethod
        def __pairwise(iterable):
            """s -> (s0,s1), (s1,s2), (s2, s3), ..."""
            a, b = tee(iterable)
            next(b, None)
            return zip(a, b)

    ReactionSupported.__qualname__ = f'ReactionSupported{_class.__name__}'
    ReactionSupported.__module__ = module or _class.__module__
    return ReactionSupported


__all__ = ['CIMtoolsTransformerMixin', 'reaction_support']
