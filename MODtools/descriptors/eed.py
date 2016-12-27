#!/usr/bin/env python3.4
# -*- coding: utf-8 -*-
#
# Copyright 2016 Ramil Nugmanov <stsouko@live.ru>
# This file is part of MODtools.
#
# MODtools is free software; you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Affero General Public License for more details.
#
#  You should have received a copy of the GNU Affero General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#
import pandas as pd
from io import StringIO
from subprocess import Popen, PIPE
from CGRtools.files.SDFrw import SDFwrite
from .basegenerator import BaseGenerator
from ..config import EED
from ..structprepare import Pharmacophoreatommarker, StandardizeDragos, CGRatommarker


class Eed(BaseGenerator):
    def __init__(self, workpath='.', s_option=None, marker_rules=None, standardize=None, cgr_reverse=False,
                 cgr_marker=None, cgr_marker_prepare=None, cgr_marker_postprocess=None, cgr_stereo=False,
                 is_reaction=False):

        if is_reaction and not cgr_marker:
            raise Exception('only cgr marker can work with reactions')

        BaseGenerator.__init__(self, workpath=workpath, s_option=s_option)

        self.__phm_marker = Pharmacophoreatommarker(marker_rules, workpath) if marker_rules else None

        self.__cgr_marker = CGRatommarker(cgr_marker, prepare=cgr_marker_prepare,
                                          postprocess=cgr_marker_postprocess,
                                          stereo=cgr_stereo, reverse=cgr_reverse) if cgr_marker else None

        self.__dragos_std = StandardizeDragos(standardize) if standardize is not None and not is_reaction else None

        self.markers = self.__cgr_marker.getcount() if cgr_marker else \
            self.__phm_marker.getcount() if marker_rules else None
        self.__workfiles = self.markers or 1

    def setworkpath(self, workpath):
        BaseGenerator.setworkpath(self, workpath)
        if self.__phm_marker:
            self.__phm_marker.setworkpath(workpath)

    def prepare(self, structures, **kwargs):
        if self.__dragos_std:
            structures = self.__dragos_std.get(structures)

        if not structures:
            return False

        if self.__cgr_marker:
            structures = self.__cgr_marker.get(structures)
        elif self.__phm_marker:
            structures = self.__phm_marker.get(structures)

        if not structures:
            return False

        workfiles = [StringIO() for _ in range(self.__workfiles)]
        writers = [SDFwrite(x, mark_to_map=True) for x in workfiles]

        prop, doubles, used_str = self.write_prepared(structures, writers)

        tX, tD = [], []
        for n, workfile in enumerate(workfiles):
            p = Popen([EED], stdout=PIPE, stdin=PIPE)
            res = p.communicate(input=workfile.getvalue().encode())[0].decode()

            if p.returncode != 0:
                return False

            X, D = self.__parseeedoutput(res)
            tX.append(X)
            tD.append(D)

        return tX, prop, tD, doubles, used_str

    @staticmethod
    def __parseeedoutput(output):
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

        x = pd.DataFrame(vector, columns=range(1, 705)).fillna(0)
        x.columns = ['eed.%d' % x for x in range(1, 705)]
        return x, pd.Series(ad)
