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
from pandas import DataFrame, Series
from subprocess import Popen, PIPE
from .basegenerator import BaseGenerator
from ..config import EED
from ..preparers.markers import PharmacophoreAtomMarker, CGRatomMarker
from ..preparers.standardizers import StandardizeDragos


class Eed(BaseGenerator):
    def __init__(self, workpath='.', s_option=None, marker_rules=None, standardize=False,
                 cgr_reverse=False, cgr_marker=None, cgr_marker_preprocess=None, cgr_marker_postprocess=None,
                 cgr_isotope=False, cgr_element=True, cgr_stereo=False, cgr_b_templates=None, cgr_m_templates=None,
                 cgr_extralabels=False, is_reaction=False):
        """
        EED wrapper.
        :param workpath: path for temp files.
        :param s_option: modeling attribute
        :param marker_rules: Dragos atom marker procedure. For molecules only. string with chemaxon Pmapper rules xml.
        :param standardize: Dragos standardization procedure. For molecules only.
          string with chemaxon standardizer rules (xml or ..- separated params)
          if False - skipped. if None - use predefined rules.
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
            raise Exception('for CGR marker is_reaction should be True')

        BaseGenerator.__init__(self, s_option=s_option)

        if marker_rules:
            self.__marker = PharmacophoreAtomMarker(marker_rules, workpath)
        elif cgr_marker:
            self.__marker = CGRatomMarker(cgr_marker, preprocess=cgr_marker_preprocess,
                                          postprocess=cgr_marker_postprocess, extralabels=cgr_extralabels,
                                          isotope=cgr_isotope, element=cgr_element, stereo=cgr_stereo,
                                          b_templates=cgr_b_templates, m_templates=cgr_m_templates, reverse=cgr_reverse)

        self.markers = (self.__marker.get_count() if self.__marker is not None else None)
        self.__workfiles = self.markers or 1

        if standardize or standardize is None:
            self.__dragos_std = StandardizeDragos(rules=standardize)

        locs = locals()
        self.__config = dict((x, locs[x]) for x in self.__optional_configs if locs[x])

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

        workfiles = [StringIO() for _ in range(self.__workfiles)]
        writers = [SDFwrite(x, mark_to_map=True) for x in workfiles]

        prop, doubles, used_str = self.write_prepared(structures, writers)

        tx, td = [], []
        for n, workfile in enumerate(workfiles):
            p = Popen([EED], stdout=PIPE, stdin=PIPE)
            res = p.communicate(input=workfile.getvalue().encode())[0].decode()

            if p.returncode != 0:
                return False

            x, d = self.__parse_eed_output(res)
            tx.append(x)
            td.append(d)

        return tx, prop, td, doubles, used_str

    @staticmethod
    def __parse_eed_output(output):
        vector, ad = [], []
        for frag in StringIO(output):
            _, *x = frag.split()
            ad.append(True)
            tmp = {}  # X vector
            for i in x:
                k, v = i.split(':')
                tmp[int(k)] = float(v.replace(',', '.'))
            if len(tmp) <= 2:
                ad[-1] = False
            vector.append(tmp)

        x = DataFrame(vector, columns=range(1, 705)).fillna(0)
        x.columns = ['eed.%d' % x for x in range(1, 705)]
        return x, Series(ad)
