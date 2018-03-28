# -*- coding: utf-8 -*-
#
#  Copyright 2015-2018 Ramil Nugmanov <stsouko@live.ru>
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
from CGRtools.containers import MoleculeContainer
from CGRtools.core import CGRcore
from operator import itemgetter
from pathlib import Path
from .chemaxon import StandardizeChemAxon
from ..common import iter2array


class StandardizeHorvat(StandardizeChemAxon):
    def __init__(self, rules=None, unwanted=None, min_ratio=2, max_ion_size=5, min_main_size=6, max_main_size=101):
        self.unwanted = self.__load_unwanted() if unwanted is None else set(unwanted)
        self.min_ratio = min_ratio
        self.max_ion_size = max_ion_size
        self.min_main_size = min_main_size
        self.max_main_size = max_main_size
        super().__init__(rules or self.__load_rules())

    def set_params(self, unwanted, **params):
        return super().set_params(unwanted=set(unwanted), **params)

    def transform(self, x):
        """
        Standardize Molecules by Dragos Horvat workflow

        step 1 (by default): dearomatize & dealkalinize, neutralize all species,
        except for FOUR-LEGGED NITROGEN, which has to be positive for else chemically incorrect
        Automatically represent N-oxides, incl. nitros, as N+-O-.
        generate major tautomer & aromatize.

        step 2: check for bizzare salts or mixtures. strip mixtures

        :param x: {array-like}, shape [n_samples] of MoleculeContainers
        :return: array of MoleculeContainers
        """

        res = []
        for s in super().transform(x):
            if s is None:
                res.append(None)
            else:
                species = sorted(((len([None for _, e in x.nodes(data='element') if e != 'H']), x)
                                  for x in CGRcore.split(s)), key=itemgetter(0))
                if species[-1][0] <= self.max_main_size and \
                        (len(species) == 1 or (species[-1][0] / species[-2][0] >= self.min_ratio and
                                               species[-2][0] <= self.max_ion_size and
                                               species[-1][0] >= self.min_main_size)) \
                        and not self.unwanted.intersection(species[-1][1]):
                    res.append(species[-1][1])
                else:
                    res.append(None)
        return iter2array(res, allow_none=True)

    @staticmethod
    def __load_rules():
        with (Path(__file__).parent / 'horvat.xml').open() as f:
            out = f.read().strip()
        return out

    @staticmethod
    def __load_unwanted():
        with (Path(__file__).parent / 'horvat.unwanted').open() as f:
            out = set(f.read().split())
        return out

    _dtype = MoleculeContainer
