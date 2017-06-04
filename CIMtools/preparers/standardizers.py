# -*- coding: utf-8 -*-
#
#  Copyright 2015-2017 Ramil Nugmanov <stsouko@live.ru>
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
from CGRtools.files.SDFrw import SDFread, SDFwrite
from io import StringIO
from json import loads
from networkx import connected_component_subgraphs
from operator import itemgetter
from os.path import join, dirname
from subprocess import Popen, PIPE, STDOUT
from . import chemax_post
from ..config import STANDARDIZER


class StandardizeDragos(object):
    def __init__(self, rules=None, unwanted=None, min_ratio=2, max_ion_size=5, min_main_size=6, max_main_size=101):
        self.__std_rules = rules or self.__load_rules()
        self.__unwanted = self.__load_unwanted() if unwanted is None else set(unwanted)
        self.__min_ratio = min_ratio
        self.__max_ion_size = max_ion_size
        self.__min_main_size = min_main_size
        self.__max_main_size = max_main_size

    def pickle(self):
        return dict(rules=self.__std_rules, unwanted=list(self.__unwanted), min_ratio=self.__min_ratio,
                    max_ion_size=self.__max_ion_size, min_main_size=self.__min_main_size,
                    max_main_size=self.__max_main_size)

    @classmethod
    def unpickle(cls, config):
        args = {'rules', 'unwanted', 'min_ratio', 'max_ion_size', 'min_main_size', 'max_main_size'}
        if args.difference(config):
            raise Exception('Invalid config')
        return StandardizeDragos(**{k: v for k, v in config.items() if k in args})

    @staticmethod
    def __load_rules():
        with open(join(dirname(__file__), "standardrules_dragos.xml")) as f:
            out = f.read().strip()
        return out

    @staticmethod
    def __load_unwanted():
        with open(join(dirname(__file__), "unwanted.elem")) as f:
            out = set(f.read().split())
        return out

    def __processor_m(self, structure):
        p = Popen([STANDARDIZER, '-c', self.__std_rules, '-f', 'SDF'], stdout=PIPE, stdin=PIPE, stderr=STDOUT)
        with StringIO() as f:
            tmp = SDFwrite(f)
            for x in structure:
                tmp.write(x)

            res = p.communicate(input=f.getvalue().encode())[0].decode()

        if p.returncode == 0:
            return SDFread(res).read()
        return False

    def __processor_s(self, structure):
        with StringIO() as f:
            SDFwrite(f).write(structure)
            data = dict(structure=f.getvalue(), parameters="smiles",
                        filterChain=[dict(filter="standardizer",
                                          parameters=dict(standardizerDefinition=self.__std_rules))])

        res = chemax_post('calculate/molExport', data)
        if res:
            res = loads(res)
            if 'isReaction' not in res:
                return SDFread(StringIO(res['structure'])).read()
        return False

    def get(self, structure):
        """
        step 1. canonical smiles, dearomatized & dealkalinized
        neutralize all species, except for FOUR-LEGGED NITROGEN, which has to be positive for else chemically incorrect
        Automatically represent N-oxides, incl. nitros, as N+-O-.
        generate major tautomer & aromatize
        """

        structure = self.__processor_m(structure) if isinstance(structure, list) else self.__processor_s(structure)
        if structure:
            """
            step 2. check for bizzare salts or mixtures
            strip mixtures
            """
            output = []
            for s in structure:
                species = sorted(((len([n for n, d in x.nodes(data=True) if d['element'] != 'H']), x) for x in
                                  connected_component_subgraphs(s)), key=itemgetter(0))
                if species[-1][0] <= self.__max_main_size \
                        and (len(species) == 1 or
                             (species[-1][0] / species[-2][0] >= self.__min_ratio and
                              species[-2][0] <= self.__max_ion_size and
                              species[-1][0] >= self.__min_main_size)) \
                        and not self.__unwanted.intersection(species[-1][1]):
                    output.append(species[-1][1])
                else:
                    return False

            return output
        return False
