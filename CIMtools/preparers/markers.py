# -*- coding: utf-8 -*-
#
#  Copyright 2017 Ramil Nugmanov <stsouko@live.ru>
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
from CGRtools.preparer import CGRcombo
from CGRtools.reactor import CGRreactor
from CGRtools.files.RDFrw import RDFread, RDFwrite
from CGRtools.files.SDFrw import SDFwrite
from io import StringIO
from itertools import product
from json import loads
from networkx import union_all
from os import close, remove
from subprocess import Popen, PIPE, STDOUT
from tempfile import mkstemp
from xml.etree import ElementTree
from . import chemax_post
from ..config import PMAPPER, STANDARDIZER


def remove_namespace(doc, namespace):
    """Remove namespace in the passed document in place."""
    ns = u'{%s}' % namespace
    nsl = len(ns)
    for elem in doc.getiterator():
        if elem.tag.startswith(ns):
            elem.tag = elem.tag[nsl:]
    return doc


class PharmacophoreAtomMarker(object):
    def __init__(self, marker_rules, workpath='.'):
        self.__marker_rules = marker_rules
        self.__markers = list(set(x.get('Symbol') for x in remove_namespace(ElementTree.fromstring(marker_rules),
                                                                            'http://www.chemaxon.com').iter('AtomSet')))
        self.set_work_path(workpath)

    def pickle(self):
        return dict(marker_rules=self.__marker_rules)

    @classmethod
    def unpickle(cls, config):
        if {'marker_rules'}.difference(config):
            raise Exception('Invalid config')
        return cls(config['marker_rules'])

    def set_work_path(self, workpath):
        self.delete_work_path()

        fd, self.__config = mkstemp(prefix='clr_', suffix='.xml', dir=workpath)
        with open(self.__config, 'w') as f:
            f.write(self.__marker_rules)
        close(fd)

    def delete_work_path(self):
        if self.__config is not None:
            remove(self.__config)
            self.__config = None

    def get_count(self):
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
                found = [[] for _ in range(self.get_count())]
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

    __config = None


class CGRatomMarker(object):
    def __init__(self, patterns, preprocess=None, postprocess=None, reverse=False, b_templates=None, m_templates=None,
                 extralabels=False, isotope=False, element=True, stereo=False):
        self.__cgr = CGRcombo(extralabels=extralabels, isotope=isotope, element=element, stereo=stereo,
                              b_templates=b_templates, m_templates=m_templates)
        self.__react = CGRreactor(stereo=stereo, extralabels=extralabels, isotope=isotope, element=element)
        self.__init_common(patterns, preprocess, postprocess, reverse)

    def _init_unpickle(self, patterns, preprocess, postprocess, reverse, **config):
        self.__cgr = CGRcombo.unpickle(dict(cgr_type='0', **config))
        self.__react = CGRreactor.unpickle(config)
        self.__init_common([ReactionContainer.unpickle(x) for x in patterns], preprocess, postprocess, reverse)

    def __init_common(self, patterns, preprocess, postprocess, reverse):
        self.__std_prerules = preprocess
        self.__std_postrules = postprocess
        self.__markers = len(patterns[0]['products'][0])
        self.__reverse = reverse
        self.__templates = patterns
        self.__patterns = self.__react.get_template_searcher(self.__react.get_templates(patterns))

    def pickle(self):
        config = self.__cgr.pickle()
        config.pop('cgr_type')
        config.update(postprocess=self.__std_postrules, preprocess=self.__std_prerules,
                      reverse=self.__reverse, patterns=[x.pickle() for x in self.__templates])
        return config

    @classmethod
    def unpickle(cls, config):
        args = {'postprocess', 'preprocess', 'reverse', 'patterns'}
        if args.difference(config):
            raise Exception('Invalid config')
        obj = cls.__new__(cls)
        obj._init_unpickle(**config)
        return obj

    def get_count(self):
        return self.__markers

    @staticmethod
    def __processor_s(structure, rules, remap=True):
        with StringIO() as f:
            RDFwrite(f).write(structure)
            res = chemax_post('calculate/molExport',
                              dict(structure=f.getvalue(), parameters="rdf",
                                   filterChain=[dict(filter="standardizer",
                                                     parameters=dict(standardizerDefinition=rules))]))
        if res:
            res = loads(res)
            if 'isReaction' in res:
                return RDFread(StringIO(res['structure']), remap=remap).read()
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
                return RDFread(StringIO(res), remap=remap).read()
        return False

    def get(self, structure):
        if self.__std_prerules:
            structure = (self.__processor_m(structure, self.__std_prerules) if isinstance(structure, list) else
                         self.__processor_s(structure, self.__std_prerules))
            if not structure:
                return False

        markslist = []
        gs = [self.__cgr.getCGR(x) for x in (structure if isinstance(structure, list) else [structure])]
        for g in gs:
            # list of list of tuples(atom, mark) of matched centers
            marks = [[[x, y['mark']] for x, y in match['products'].nodes(data=True)] for match in self.__patterns(g)]
            markslist.append(marks)

        if self.__std_postrules:
            structure = (self.__processor_m(structure, self.__std_postrules, remap=False)
                         if isinstance(structure, list) else
                         self.__processor_s(structure, self.__std_postrules, remap=False))

            if not structure:
                return False

        output = []
        for s, marks in zip((structure if isinstance(structure, list) else [structure]), markslist):
            ss = union_all(s.products if self.__reverse else s.substrats)
            ss.meta.update(s.meta)
            ps = []
            for x in (s.products if self.__reverse else s.substrats):
                tmp = x.copy()
                tmp.meta.update(s.meta)
                ps.append(tmp)

            result = []
            for match in marks:
                tmp = []
                for atom, a_mark in match:
                    ssc = next(x for x in ps if atom in x).copy()
                    ssc.node[atom]['mark'] = '1'
                    tmp.append([a_mark, atom, ssc])

                result.append([(x, y) for _, x, y in sorted(tmp)])

            output.append(result if result else [[(None, ss)] * self.__markers])
        return output
