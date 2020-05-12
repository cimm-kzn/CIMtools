# -*- coding: utf-8 -*-
#
#  Copyright 2020 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from CGRtools import CGRContainer, MoleculeContainer
from collections.abc import Sequence
from numpy import zeros
from ..base import CIMtoolsTransformerMixin


class SlicedTuple(Sequence):
    def __init__(self, data):
        self.data = data

    def __getitem__(self, i):
        if isinstance(i, slice):
            return SlicedTuple(tuple(x[i] for x in self.data))
        else:
            return tuple(x[i] for x in self.data)

    def __len__(self):
        return len(self.data[0])

    def __repr__(self):
        return repr(self.data)


class GraphToMatrix(CIMtoolsTransformerMixin):
    def transform(self, x):
        x = super().transform(x)

        size = max(len(g) for g in x)
        atoms = zeros((len(x), size, self._atom_vector_size), dtype=int)
        bonds = zeros((len(x), size, size), dtype=int)

        for i, g in enumerate(x):
            aam = {}
            for j, (n, a) in enumerate(g.atoms()):
                aam[n] = j
                atoms[i, j, :] = self._atom_vector(a)
            for n, m, b in g.bonds():
                bonds[i, aam[n], aam[m]] = bonds[i, aam[m], aam[n]] = int(b)
        return SlicedTuple((atoms, bonds))


class MoleculesToMatrix(GraphToMatrix):
    def __init__(self, charge=True, is_radical=False, isotope=False, hybridization=False, neighbors=False,
                 implicit_hydrogens=False, total_hydrogens=False, in_ring=False):
        self.charge = charge
        self.is_radical = is_radical
        self.isotope = isotope
        self.hybridization = hybridization
        self.neighbors = neighbors
        self.implicit_hydrogens = implicit_hydrogens
        self.total_hydrogens = total_hydrogens
        self.in_ring = in_ring

    def _atom_vector(self, atom):
        vector = [atom.atomic_number]
        if self.charge:
            vector.append(atom.charge)
        if self.is_radical:
            vector.append(int(atom.is_radical))
        if self.isotope:
            vector.append(atom.isotope or 0)
        if self.hybridization:
            vector.append(atom.hybridization)
        if self.neighbors:
            vector.append(atom.neighbors)
        if self.implicit_hydrogens:
            vector.append(atom.implicit_hydrogens)
        if self.total_hydrogens:
            vector.append(atom.total_hydrogens)
        if self.in_ring:
            vector.append(int(atom.in_ring))
        return vector

    @property
    def _atom_vector_size(self):
        size = 1
        if self.charge:
            size += 1
        if self.is_radical:
            size += 1
        if self.isotope:
            size += 1
        if self.hybridization:
            size += 1
        if self.neighbors:
            size += 1
        if self.implicit_hydrogens:
            size += 1
        if self.total_hydrogens:
            size += 1
        if self.in_ring:
            size += 1
        return size

    _dtype = MoleculeContainer


class CGRToMatrix(GraphToMatrix):
    def __init__(self, charge=True, is_radical=False, isotope=False, hybridization=False, neighbors=False,
                 in_ring=False):
        self.charge = charge
        self.is_radical = is_radical
        self.isotope = isotope
        self.hybridization = hybridization
        self.neighbors = neighbors
        self.in_ring = in_ring

    def _atom_vector(self, atom):
        vector = [atom.atomic_number]
        if self.charge:
            vector.append(atom.charge)
            vector.append(atom.p_charge)
        if self.is_radical:
            vector.append(int(atom.is_radical))
            vector.append(int(atom.p_is_radical))
        if self.isotope:
            vector.append(atom.isotope or 0)
        if self.hybridization:
            vector.append(atom.hybridization)
            vector.append(atom.p_hybridization)
        if self.neighbors:
            vector.append(atom.neighbors)
            vector.append(atom.p_neighbors)
        if self.in_ring:
            vector.append(int(atom.in_ring))
        return vector

    @property
    def _atom_vector_size(self):
        size = 1
        if self.charge:
            size += 2
        if self.is_radical:
            size += 2
        if self.isotope:
            size += 1
        if self.hybridization:
            size += 2
        if self.neighbors:
            size += 2
        if self.in_ring:
            size += 1
        return size

    _dtype = CGRContainer


__all__ = ['MoleculesToMatrix', 'CGRToMatrix']
