# -*- coding: utf-8 -*-
#
# Copyright 2015, 2016 Ramil Nugmanov <stsouko@live.ru>
# This file is part of MODtools.
#
# MODtools is free software; you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Affero General Public License for more details.
#
#  You should have received a copy of the GNU Affero General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#
import json
import networkx as nx
import os
import xml.etree.ElementTree as ET
from io import StringIO
from itertools import product
from operator import itemgetter
from subprocess import Popen, PIPE, STDOUT, call
from CGRtools.CGRpreparer import CGRbalanser, CGRcombo
from CGRtools.CGRreactor import CGRreactor
from CGRtools.files.RDFrw import RDFread, RDFwrite
from CGRtools.files.SDFrw import SDFread, SDFwrite
from .config import PMAPPER, STANDARDIZER, COLOR
from .utils import chemaxpost
from . import remove_namespace


class StandardizeDragos(object):
    def __init__(self, rules):
        self.__stdrules = self.__dumprules(rules)
        self.__unwanted = self.__loadunwanted()
        self.__minratio = 2
        self.__maxionsize = 5
        self.__minmainsize = 6
        self.__maxmainsize = 101

    def __dumprules(self, rules):
        with rules or open(os.path.join(os.path.dirname(__file__), "standardrules_dragos.rules")) as f:
            ruless = f.read()
        return ruless

    def __loadunwanted(self):
        return set(open(os.path.join(os.path.dirname(__file__), "unwanted.elem")).read().split())

    def __processor_m(self, structure):
        p = Popen([STANDARDIZER, '-c', self.__stdrules, '-f', 'SDF'], stdout=PIPE, stdin=PIPE, stderr=STDOUT)
        with StringIO() as f:
            tmp = SDFwrite(f)
            for x in structure:
                tmp.write(x)

            res = p.communicate(input=f.getvalue().encode())[0].decode()

        if p.returncode == 0:
            with StringIO(res) as f:
                return list(SDFread(f).read())
        return False

    def __processor_s(self, structure):
        with StringIO() as f:
            SDFwrite(f).write(structure)
            data = {"structure": f.getvalue(), "parameters": "smiles",
                    "filterChain": [{"filter": "standardizer",
                                     "parameters": {"standardizerDefinition": self.__stdrules}}]}

        res = chemaxpost('calculate/molExport', data)
        if res:
            res = json.loads(res)
            if 'isReaction' not in res:
                with StringIO(res['structure']) as f:
                    return list(SDFread(f).read())
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
                                  nx.connected_component_subgraphs(s)), key=itemgetter(0))
                if species[-1][0] <= self.__maxmainsize \
                        and (len(species) == 1 or
                             (species[-1][0] / species[-2][0] >= self.__minratio and
                              species[-2][0] <= self.__maxionsize and
                              species[-1][0] >= self.__minmainsize)) \
                        and not self.__unwanted.intersection(species[-1][1]):
                    output.append(species[-1][1])
                else:
                    return False

            return output
        return False


class Pharmacophoreatommarker(object):
    def __init__(self, markerrule, workpath):
        self.__markerrule, self.__markers = self.__dumprules(markerrule)
        self.setworkpath(workpath)

    @staticmethod
    def __dumprules(rules):
        with rules as f:
            _rules = f.read()
        marks = list(set(x.get('Symbol') for x in remove_namespace(ET.fromstring(_rules),
                                                                   'http://www.chemaxon.com').iter('AtomSet')))
        return _rules, marks

    def setworkpath(self, workpath):
        self.__config = os.path.join(workpath, 'iam')
        with open(self.__config, 'w') as f:
            f.write(self.__markerrule)

    def getcount(self):
        return len(self.__markers)

    def get(self, structure):
        """
        marks atoms in 7th col of sdf.
        if molecule has atom mapping - will be used mapping.
        :type structure: nx.Graph or list(nx.Graph)
        """
        p = Popen([PMAPPER, '-c', self.__config], stdout=PIPE, stdin=PIPE, stderr=STDOUT)
        with StringIO() as f:
            tmp = SDFwrite(f)
            for x in (structure if isinstance(structure, list) else [structure]):
                tmp.write(x)

            marks = p.communicate(input=f.getvalue().encode())[0].decode().split()

        if p.returncode == 0:
            output = []
            for s, mark in zip((structure if isinstance(structure, list) else [structure]), marks):
                found = [[] for _ in range(self.getcount())]
                for n, m in zip(s.nodes(), mark.split(';')):
                    if m:
                        tmp = s.copy()
                        tmp.node[n]['mark'] = '1'
                        found[self.__markers.index(m)].append((n, tmp))

                for x in found:
                    if not x:
                        x.append((None, s))

                output.append([list(x) for x in product(*found)])
            return output
        return False


class CGRatommarker(CGRcombo, CGRreactor):
    def __init__(self, patterns, prepare=None, postprocess=None, reverse=False,
                 b_templates=None, m_templates=None, speed=False, extralabels=False, isotop=False, element=True, deep=0,
                 stereo=False):

        CGRreactor.__init__(self, stereo=stereo, hyb=extralabels, neighbors=extralabels, isotop=isotop, element=element,
                            deep=deep)

        CGRcombo.__init__(self, cgr_type='0',
                          extralabels=extralabels, isotop=isotop, element=element, deep=deep, stereo=stereo,
                          b_templates=b_templates, m_templates=m_templates, speed=speed)

        self.__stdprerules = self.__dumprules(prepare)
        self.__stdpostrules = self.__dumprules(postprocess)
        self.__templates, self.__marks = self.__dumppatterns(patterns)
        self.__reverse = reverse

    @staticmethod
    def __dumprules(rules):
        if rules:
            with rules as f:
                rules = f.read().rstrip()
        return rules

    @staticmethod
    def __dumppatterns(patterns):
        with patterns as f:
            templates = CGRbalanser.get_templates(f)
            marks = len(templates[0]['products'])
        return templates, marks

    def getcount(self):
        return self.__marks

    @staticmethod
    def __processor_s(structure, rules, remap=True):
        with StringIO() as f:
            RDFwrite(f).write(structure)
            res = chemaxpost('calculate/molExport',
                             {"structure": f.getvalue(), "parameters": "rdf",
                              "filterChain": [{"filter": "standardizer",
                                               "parameters": {"standardizerDefinition": rules}}]})
        if res:
            res = json.loads(res)
            if 'isReaction' in res:
                return list(RDFread(StringIO(res['structure'])).read(remap=remap))
        return False

    @staticmethod
    def __processor_m(structure, rules, remap=True):
        p = Popen([STANDARDIZER, '-c', rules, '-f', 'rdf'], stdout=PIPE, stdin=PIPE, stderr=STDOUT)
        with StringIO() as f:
            tmp = RDFwrite(f)
            for x in structure:
                tmp.write(x)

            res = p.communicate(input=f.getvalue().encode())[0].decode()
            if p.returncode == 0:
                return list(RDFread(StringIO(res)).read(remap=remap))
        return False

    def get(self, structure):
        if self.__stdprerules:
            structure = self.__processor_m(structure, self.__stdprerules) if isinstance(structure, list) \
                else self.__processor_s(structure, self.__stdprerules)
            if not structure:
                return False

        _patterns = self.searchtemplate(self.__templates, speed=False)  # ad_hoc for pickle

        markslist = []
        gs = [self.getCGR(x) for x in (structure if isinstance(structure, list) else [structure])]
        for g in gs:
            # list of list of tuples(atom, mark) of matched centers
            marks = [[[x, y['mark']] for x, y in match['products'].nodes(data=True)] for match in _patterns(g)]
            markslist.append(marks)

        if self.__stdpostrules:
            structure = self.__processor_m(structure, self.__stdpostrules, remap=False) if isinstance(structure, list) \
                else self.__processor_s(structure, self.__stdpostrules, remap=False)

            if not structure:
                return False

        output = []
        for s, marks in zip((structure if isinstance(structure, list) else [structure]), markslist):
            ss = nx.union_all(s['products'] if self.__reverse else s['substrats'])
            ss.graph['meta'] = s['meta'].copy()
            ps = []
            for x in (s['products'] if self.__reverse else s['substrats']):
                tmp = x.copy()
                tmp.graph['meta'] = s['meta'].copy()
                ps.append(tmp)

            result = []
            for match in marks:
                tmp = []
                for atom, a_mark in match:
                    ssc = next(x for x in ps if atom in x).copy()
                    ssc.node[atom]['mark'] = '1'
                    tmp.append([a_mark, atom, ssc])

                result.append([(x, y) for _, x, y in sorted(tmp)])

            output.append(result if result else [[(None, ss)] * self.getcount()])
        return output


class Colorize(object):
    def __init__(self, standardize, workpath):
        self.__standardize = self.__dumprules(standardize)
        self.setworkpath(workpath)

    @staticmethod
    def __dumprules(rules):
        with rules or open(os.path.join(os.path.dirname(__file__), "standardrules_dragos.rules")) as f:
            rules = f.read()
        return rules

    def setworkpath(self, workpath):
        self.__input_file = os.path.join(workpath, 'colorin.sdf')
        self.__out_file = os.path.join(workpath, 'colorout.sdf')
        self.__std_file = os.path.join(workpath, 'colorstd.xml')
        with open(self.__std_file, 'w') as f:
            f.write(self.__standardize)

    def get(self, structure):
        if os.path.exists(self.__out_file):
            os.remove(self.__out_file)
        with open(self.__input_file, 'w') as f:
            out = SDFwrite(f)
            for i in (structure if isinstance(structure, list) else [structure]):
                out.write(i)

        if call([COLOR, self.__input_file, self.__out_file, self.__std_file]) == 0:
            with open(self.__out_file) as f:
                res = list(SDFread(f).read(remap=False))
                if res:
                    return res
        return False
