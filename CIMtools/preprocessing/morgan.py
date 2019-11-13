# -*- coding: utf-8 -*-
#
#  Copyright 2019 Ramil Nugmanov <stsouko@live.ru>
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
from CGRtools.containers import MoleculeContainer, CGRContainer, QueryContainer, QueryCGRContainer
from collections import  defaultdict
from numpy import zeros, bool8
from hashlib import md5
from sklearn.base import BaseEstimator, TransformerMixin
from ..utils import iter2array


class MorganFingerprint(BaseEstimator, TransformerMixin):
    def __init__(self, fingerprint_size=12, deep=1):
        self.fingerprint_size = fingerprint_size
        self.deep = deep

    def fit(self, x, y=None):
        return self

    def transform(self, x):
        return x

    def _transform_bitset(self, x):
        pass

    def _transform_hashes(self, x):
        x = iter2array(x, dtype=(MoleculeContainer, CGRContainer, QueryContainer, QueryCGRContainer))
        results = []
        for mol in x:
            adj = defaultdict(dict)
            for n, m, b in mol.bonds():
                adj[n][m] = adj[m][n] = int(b)
            nodes = {n: int(a) for n, a in mol.atoms()}
            hash_set = set(nodes.values())
            for _ in range(self.deep):
                nodes = {n: hash((a, tuple(sorted((b, nodes[m]) for m, b in adj[n])))) for n, a in nodes.items()}
                hash_set.update(nodes.values())
            results.append(hash_set)
        return results