#!/usr/bin/env python3.4
# -*- coding: utf-8 -*-
#
#  Copyright 2016 Ramil Nugmanov <stsouko@live.ru>
#  This file is part of MODtools.
#
#  MODtools
#  is free software; you can redistribute it and/or modify
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
import pandas as pd
from io import StringIO
from subprocess import Popen, PIPE, STDOUT
from CGRtools.files.RDFrw import RDFread
from CGRtools.files.SDFrw import SDFread, SDFwrite
from ..config import CXCALC
from .descriptoragregator import Propertyextractor
from ..structprepare import Pharmacophoreatommarker, StandardizeDragos, CGRatommarker


class Pkab(Propertyextractor):
    def __init__(self, workpath='.', s_option=None, marker_rules=None, standardize=None, acid=True, base=True,
                 cgr_reverse=False, is_reaction=False,
                 cgr_marker=None, cgr_marker_prepare=None, cgr_marker_postprocess=None, cgr_stereo=False):

        self.__is_reaction = is_reaction
        if is_reaction and not cgr_marker:
            raise Exception('only cgr marker can work with reactions')

        Propertyextractor.__init__(self, s_option)

        self.__dragos_marker = Pharmacophoreatommarker(marker_rules, workpath) if marker_rules else None

        self.__cgr_marker = CGRatommarker(cgr_marker, prepare=cgr_marker_prepare,
                                          postprocess=cgr_marker_postprocess,
                                          stereo=cgr_stereo, reverse=cgr_reverse) if cgr_marker else None

        self.__dragos_std = StandardizeDragos(standardize) \
            if standardize is not None and not self.__cgr_marker else None
        self.__workpath = workpath
        self.__reverse = cgr_reverse
        self.__acid = acid
        self.__base = base

    def setworkpath(self, workpath):
        self.__workpath = workpath
        if self.__dragos_marker:
            self.__dragos_marker.setworkpath(workpath)

    def get(self, structures, **kwargs):
        reader = RDFread(structures) if self.__is_reaction else SDFread(structures)
        data = list(reader.read())
        structures.seek(0)  # ad-hoc for rereading

        if self.__dragos_std:
            data = self.__dragos_std.get(data)

        if not data:
            return False

        if self.__cgr_marker:
            data = self.__cgr_marker.get(data)

        elif self.__dragos_marker:
            data = self.__dragos_marker.get(data)

        if not data:
            return False

        prop = []
        doubles = []

        p = Popen([CXCALC] + 'pka -x 50 -i -50 -a 8 -b 8 -P dynamic -m micro'.split(),
                  stdout=PIPE, stdin=PIPE, stderr=STDOUT)

        with StringIO() as f:
            writer = SDFwrite(f)
            for s_numb, s in enumerate(data):
                if self.__cgr_marker:
                    meta = s[0][0][1].graph['meta']
                    for d in s:
                        tmp = [s_numb]
                        for x, y in d:
                            writer.write(y)
                            tmp.append(x)
                        prop.append(self.get_property(meta, marks=tmp[1:]))
                        doubles.extend([tmp] * len(d))
                elif self.__dragos_marker:
                    writer.write(s[0][0][1])
                    prop.extend([self.get_property(d[0][1].graph['meta'], marks=[x[0] for x in d]) for d in s])
                    doubles.append([([s_numb] + [x[0] for x in d]) for d in s])
                else:
                    writer.write(s)
                    prop.append(self.get_property(s.graph['meta']))
                    doubles.append(s_numb)

            res = p.communicate(input=f.getvalue().encode())[0].decode()

        if p.returncode == 0:
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
            if self.__cgr_marker or self.__dragos_marker:
                old_k = None
                for lk, v in zip(doubles, pk):
                    for k in (lk if self.__dragos_marker else [lk]):
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

            X = pd.DataFrame([x[1] for x in new_doubles],
                             columns=(['pka'] if self.__acid else []) + (['pkb'] if self.__base else []))

            res = dict(X=X, AD=-X.isnull().any(axis=1), Y=pd.Series(prop, name='Property'))

            if self.__cgr_marker or self.__dragos_marker:
                i = pd.MultiIndex.from_tuples([x[0] for x in new_doubles], names=['structure', 'c.0', 'c.1'])
            else:
                i = pd.Index(doubles, name='structure')

            res['X'].index = res['AD'].index = res['Y'].index = i
            return res

        return False
