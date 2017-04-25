# -*- coding: utf-8 -*-
#
#  Copyright 2016, 2017 Ramil Nugmanov <stsouko@live.ru>
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
from CGRtools.files.SDFrw import SDFwrite
from io import StringIO
from pandas import DataFrame
from subprocess import Popen, PIPE, STDOUT
from .basegenerator import BaseGenerator
from ..config import CXCALC
from ..preparers.markers import PharmacophoreAtomMarker, CGRatomMarker
from ..preparers.standardizers import StandardizeDragos


class Pkab(BaseGenerator):
    def __init__(self, workpath='.', s_option=None, marker_rules=None, standardize=False, acid=True, base=True,
                 cgr_reverse=False, cgr_marker=None, cgr_marker_preprocess=None, cgr_marker_postprocess=None,
                 cgr_isotope=False, cgr_element=True, cgr_stereo=False, cgr_b_templates=None, cgr_m_templates=None,
                 cgr_extralabels=False, is_reaction=False):
        """
        Chemaxon cxcalc wrapper
        :param workpath: path for temp files.
        :param s_option: modeling attribute
        :param marker_rules: Dragos atom marker procedure. For molecules only. string with chemaxon Pmapper rules xml.
        :param standardize: Dragos standardization procedure. For molecules only.
          string with chemaxon standardizer rules (xml or ..- separated params)
        :param acid: calc acidity
        :param base: calc basicity
        :param cgr_reverse: for reactions only. use product structure for calculation.
        :param cgr_marker: CGRtools.CGRreactor templates  
        :param cgr_marker_preprocess: prepare (if need) reaction structure for marker. 
          string with chemaxon standardizer rules (xml or ..- separated params)
        :param cgr_marker_postprocess: postprocess (if need) after marker [previously prepared] reaction structure  
          string with chemaxon standardizer rules (xml or ..- separated params)
        :param cgr_b_templates: list of reactioncontainers with termplates for autobalancing reactions
        :param cgr_m_templates: list of reactioncontainers with termplates for reactions mapping correction
        :param cgr_extralabels: match neighbors and hyb marks in substructure search. 
          Need for CGRatomMarker or CGRcombo balanser/remapper
        :param cgr_isotope: match isotope. see cgr_extralabels
        :param cgr_element: match stereo. see cgr_extralabels
        :param cgr_stereo: match stereo.
        :param is_reaction: True for reaction data
        """
        if is_reaction:
            if not cgr_marker:
                raise Exception('need CGR marker for work with reactions')
            if standardize or standardize is None:
                raise Exception('standardize can work only with molecules')
            if marker_rules:
                raise Exception('pharmacophore atom marker can work only with molecules')
        elif cgr_marker:
            raise Exception('for cgr marker is_reaction should be True')

        BaseGenerator.__init__(self, s_option=s_option)

        if marker_rules:
            self.__marker = PharmacophoreAtomMarker(marker_rules, workpath)
        elif cgr_marker:
            self.__marker = CGRatomMarker(cgr_marker, preprocess=cgr_marker_preprocess,
                                          postprocess=cgr_marker_postprocess, extralabels=cgr_extralabels,
                                          isotope=cgr_isotope, element=cgr_element, stereo=cgr_stereo,
                                          b_templates=cgr_b_templates, m_templates=cgr_m_templates, reverse=cgr_reverse)

        self.markers = (self.__marker.get_count() if self.__marker is not None else None)

        if standardize or standardize is None:
            self.__dragos_std = StandardizeDragos(rules=standardize)

        self.__reverse = cgr_reverse
        self.__acid = acid
        self.__base = base

        tmp = []
        if acid:
            tmp.append('pka')
        if base:
            tmp.append('pkb')
        self.__columns = tmp

        locs = locals()
        tmp = dict((x, locs[x]) for x in self.__optional_configs if locs[x])
        tmp.update((x, y) for x, y in (('acid', acid), ('base', base)) if not y)
        self.__config = tmp

    __optional_configs = ('s_option', 'marker_rules', 'standardize', 'cgr_reverse', 'cgr_marker',
                          'cgr_marker_preprocess', 'cgr_marker_postprocess', 'cgr_stereo', 'is_reaction')
    __marker = None
    __dragos_std = None

    def get_config(self):
        return self.__config

    def set_work_path(self, workpath):
        if hasattr(self.__marker, 'set_work_path'):
            self.__marker.set_work_path(workpath)

    def delete_work_path(self):
        if hasattr(self.__marker, 'delete_work_path'):
            self.__marker.delete_work_path()

    def pickle(self):
        if hasattr(self.__marker, 'pickle'):
            self.__marker.pickle()

    def prepare(self, structures, **_):
        if self.__dragos_std:
            structures = self.__dragos_std.get(structures)

        if not structures:
            return False

        if self.__marker:
            structures = self.__marker.get(structures)

        if not structures:
            return False

        prop, doubles, used_str = [], [], []

        p = Popen([CXCALC] + 'pka -x 50 -i -50 -a 8 -b 8 -P dynamic -m micro'.split(),
                  stdout=PIPE, stdin=PIPE, stderr=STDOUT)

        with StringIO() as f:
            writer = SDFwrite(f)
            for s_numb, s in enumerate(structures):
                if isinstance(self.__marker, CGRatomMarker):
                    meta = s[0][0][1].meta
                    for d in s:
                        tmp_d = [s_numb]
                        tmp_s = []  # list of graphs with marked atoms
                        for x, y in d:
                            writer.write(y)
                            tmp_d.append(x)
                            tmp_s.append(y)
                        prop.append(self.get_property(meta, marks=tmp_d[1:]))
                        doubles.extend([tmp_d] * len(d))
                        used_str.append(tmp_s)
                elif isinstance(self.__marker, PharmacophoreAtomMarker):
                    writer.write(s[0][0][1])
                    prop.extend([self.get_property(d[0][1].meta, marks=[x[0] for x in d]) for d in s])
                    doubles.append([([s_numb] + [x[0] for x in d]) for d in s])
                    used_str.extend([s[0][0][1]]*len(s))
                else:
                    writer.write(s)
                    prop.append(self.get_property(s.meta))
                    doubles.append(s_numb)
                    used_str.append(s)

            res = p.communicate(input=f.getvalue().encode())[0].decode()

        if p.returncode != 0:
            return False

        pk = []
        with StringIO(res) as f:
            f.readline()
            for i in f:
                ii = i.rstrip().split('\t')
                key = iter(ii[-1].split(','))
                val = ii[1: -1]
                pk.append([{int(next(key)): float(x.replace(',', '.')) for x in v if x}
                           for v in [val[:8], val[8:]]])

        new_doubles = []
        if self.__marker:
            old_k = None
            for lk, v in zip(doubles, pk):
                for k in (lk if isinstance(self.__marker, PharmacophoreAtomMarker) else [lk]):
                    if k == old_k:
                        if self.__base:
                            new_doubles[-1][1].append(v[1].get(k[1 if self.__reverse else 2]))
                    else:
                        old_k = k
                        new_doubles.append([k, ([v[0].get(k[2 if self.__reverse else 1])] if self.__acid else [])])

        else:
            for k, v in zip(doubles, pk):
                new_doubles.append([k, ([min(v[0].values())] if self.__acid else []) +
                                       ([max(v[1].values())] if self.__base else [])])

        x = DataFrame([x[1] for x in new_doubles], columns=self.__columns)
        ad = [x.isnull().any(axis=1) ^ True]
        i = [x[0] for x in new_doubles] if self.markers else doubles

        return [x], prop, [ad], i, used_str

    def write_prepared(self, *_):
        raise Exception('Disabled method')
