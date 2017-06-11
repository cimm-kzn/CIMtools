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


class Eed(BaseGenerator):
    def __init__(self, workpath='.', s_option=None, marker_rules=None, standardize=False,
                 cgr_reverse=False, cgr_marker=None, cgr_marker_preprocess=None, cgr_marker_postprocess=None,
                 cgr_isotope=False, cgr_element=True, cgr_stereo=False, cgr_b_templates=None, cgr_m_templates=None,
                 cgr_extralabels=False, is_reaction=False):
        """
        EED wrapper.
        """
        if is_reaction and not cgr_marker:
            raise Exception('need CGR marker for work with reactions')

        BaseGenerator.__init__(self, workpath=workpath, s_option=s_option, marker_rules=marker_rules,
                               standardize=standardize, cgr_marker=cgr_marker, cgr_isotope=cgr_isotope,
                               cgr_marker_preprocess=cgr_marker_preprocess, cgr_extralabels=cgr_extralabels,
                               cgr_marker_postprocess=cgr_marker_postprocess, cgr_element=cgr_element,
                               cgr_stereo=cgr_stereo, cgr_b_templates=cgr_b_templates, cgr_m_templates=cgr_m_templates,
                               cgr_reverse=cgr_reverse, is_reaction=is_reaction)

        self.__workfiles = self._markers_count or 1

    def _init_unpickle(self, s_option, cgr_isotope, cgr_element, cgr_stereo, cgr_marker_preprocess,
                       cgr_marker_postprocess, cgr_extralabels, cgr_reverse, is_reaction, marker_rules, cgr_marker,
                       cgr_b_templates, cgr_m_templates, standardize, **config):
        BaseGenerator._init_unpickle(self, s_option, cgr_isotope, cgr_element, cgr_stereo, cgr_marker_preprocess,
                                     cgr_marker_postprocess, cgr_extralabels, cgr_reverse, is_reaction,
                                     marker_rules, cgr_marker, cgr_b_templates, cgr_m_templates, standardize, **config)
        self.__workfiles = self._markers_count or 1

    def pickle(self):
        return super(Eed, self).pickle()

    @classmethod
    def unpickle(cls, config):
        BaseGenerator.unpickle(config)
        obj = cls.__new__(cls)
        obj._init_unpickle(**config)
        return obj

    def set_work_path(self, workpath):
        super(Eed, self).set_work_path()

    def delete_work_path(self):
        super(Eed, self).delete_work_path()

    def _prepare(self, structures, **_):
        if self._dragos_std:
            structures = self._dragos_std.get(structures)

        if not structures:
            return False

        if self._marker:
            structures = self._marker.get(structures)

        if not structures:
            return False

        workfiles = [StringIO() for _ in range(self.__workfiles)]
        writers = [SDFwrite(x, mark_to_map=True) for x in workfiles]

        prop, doubles = self._write_prepared(structures, writers)

        tx, td = [], []
        for n, workfile in enumerate(workfiles):
            p = Popen([EED], stdout=PIPE, stdin=PIPE)
            res = p.communicate(input=workfile.getvalue().encode())[0].decode()

            if p.returncode != 0:
                return False

            x, d = self.__parse_eed_output(res)
            tx.append(x)
            td.append(d)

        return tx, prop, td, doubles

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
