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
from abc import ABC, abstractmethod
from functools import reduce
from operator import and_
from pandas import Series, Index, MultiIndex, concat
from .propertyextractor import PropertyExtractor
from ..preparers.markers import PharmacophoreAtomMarker, CGRatomMarker
from ..preparers.standardizers import StandardizeDragos


class BaseGenerator(ABC, PropertyExtractor):
    def __init__(self, workpath='.', s_option=None, marker_rules=None, standardize=False, cgr_marker=None,
                 cgr_marker_preprocess=None, cgr_marker_postprocess=None,  cgr_extralabels=False, cgr_isotope=False,
                 cgr_element=True, cgr_stereo=False, cgr_b_templates=None, cgr_m_templates=None, cgr_reverse=False,
                 is_reaction=False):
        """
        :param workpath: path for temp files.
        :param s_option: modeling attribute
        :param marker_rules: Dragos atom marker procedure. For molecules only. string with chemaxon Pmapper rules xml.
        :param standardize: Dragos standardization procedure. For molecules only.
          string with chemaxon standardizer rules (xml or ..- separated params)
          if False - skipped. if None - use predefined rules.
        :param cgr_marker: CGRtools.CGRreactor templates  
        :param cgr_marker_preprocess: prepare (if need) reaction structure for marker. 
          string with chemaxon standardizer rules (xml or ..- separated params)
        :param cgr_marker_postprocess: postprocess (if need) after marker [previously prepared] reaction structure  
          string with chemaxon standardizer rules (xml or ..- separated params)
        :param cgr_extralabels: match neighbors and hyb marks in substructure search. 
          Need for CGRatomMarker or CGRcombo balanser/remapper
        
        :param cgr_isotope: match isotope. see cgr_extralabels
        :param cgr_element: match stereo. see cgr_extralabels
        :param cgr_stereo: match stereo.
        :param cgr_b_templates: list of reactioncontainers with termplates for autobalancing reactions
        :param cgr_m_templates: list of reactioncontainers with termplates for reactions mapping correction
        :param cgr_reverse: for reactions only. use product structure for calculation.
        :param is_reaction: True for reaction data
        """
        if is_reaction:
            if not cgr_marker:
                raise Exception('only cgr or cgr marker can work with reactions')
            if standardize or standardize is None:
                raise Exception('standardize can work only with molecules')
            if marker_rules:
                raise Exception('pharmacophore atom marker can work only with molecules')
        elif cgr_marker:
            raise Exception('for cgr_marker is_reaction should be True')

        PropertyExtractor.__init__(self, s_option)

        if standardize or standardize is None:
            self._dragos_std = StandardizeDragos(rules=standardize)

        if marker_rules:
            self._marker = PharmacophoreAtomMarker(marker_rules, workpath)
        elif cgr_marker:
            self._marker = CGRatomMarker(cgr_marker, preprocess=cgr_marker_preprocess,
                                         postprocess=cgr_marker_postprocess, extralabels=cgr_extralabels,
                                         isotope=cgr_isotope, element=cgr_element, stereo=cgr_stereo,
                                         b_templates=cgr_b_templates, m_templates=cgr_m_templates, reverse=cgr_reverse)

        self._markers_count = self._marker.get_count() if self._marker is not None else None
        self.__pickle = dict(s_option=s_option, cgr_isotope=cgr_isotope, cgr_element=cgr_element, cgr_stereo=cgr_stereo,
                             cgr_marker_preprocess=cgr_marker_preprocess, cgr_marker_postprocess=cgr_marker_postprocess,
                             cgr_extralabels=cgr_extralabels, cgr_reverse=cgr_reverse, is_reaction=is_reaction,
                             marker_rules=None, cgr_marker=None, cgr_b_templates=None, cgr_m_templates=None,
                             standardize=False)

    @abstractmethod
    def pickle(self):
        config = self.__pickle.copy()
        if isinstance(self._marker, PharmacophoreAtomMarker):
            config.update(self._marker.pickle())
        elif isinstance(self._marker, CGRatomMarker):
            x = self._marker.pickle()
            config.update(cgr_marker=x['patterns'], cgr_b_templates=x['b_templates'], cgr_m_templates=x['m_templates'])
        if self._dragos_std is not None:
            config['standardize'] = self._dragos_std.pickle()['rules']
        return config

    @abstractmethod
    def prepare(self, structures, **_):
        pass

    @abstractmethod
    def set_work_path(self, workpath):
        pass

    @abstractmethod
    def delete_work_path(self):
        pass

    def get(self, structures, **kwargs):
        tmp = self.prepare(structures, **kwargs)
        if not tmp:
            return False

        x, y, ad, i, s = tmp

        res = dict(X=concat(x, axis=1, keys=range(len(x))) if len(x) > 1 else x[0],
                   AD=reduce(and_, ad), Y=Series(y, name='Property'), structures=s)  # todo: prepare structures

        if self._markers_count:
            _i = MultiIndex.from_tuples(i, names=['structure'] + ['c.%d' % x for x in range(self._markers_count)])
        else:
            _i = Index(i, name='structure')

        res['X'].index = res['AD'].index = res['Y'].index = _i
        return res

    def write_prepared(self, structures, writers):
        prop = []
        doubles = []
        used_str = []
        for s_numb, s in enumerate(structures):
            if isinstance(s, list):
                meta = s[0][0][1].meta
                for d in s:  # d = ((n1, tmp1), (n2, tmp2), ...)
                    tmp_d = [s_numb]
                    tmp_s = []  # list of graphs with marked atoms
                    for w, (x, y) in zip(writers, d):
                        w.write(y)
                        tmp_d.append(x)
                        tmp_s.append(y)
                    prop.append(self.get_property(meta, marks=tmp_d[1:]))
                    doubles.append(tmp_d)
                    used_str.append(tmp_s)
            else:
                writers[0].write(s)
                prop.append(self.get_property(s.meta))
                doubles.append(s_numb)
                used_str.append(s)

        return prop, doubles, used_str

    _marker = _dragos_std = None
