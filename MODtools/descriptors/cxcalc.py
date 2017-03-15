# -*- coding: utf-8 -*-
#
#  Copyright 2016, 2017 Ramil Nugmanov <stsouko@live.ru>
#  This file is part of MODtools.
#
#  MODtools is free software; you can redistribute it and/or modify
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
from pandas import DataFrame, Series, Index, MultiIndex
from io import StringIO
from subprocess import Popen, PIPE, STDOUT
from CGRtools.files.SDFrw import SDFwrite
from ..config import CXCALC
from .descriptoragregator import PropertyExtractor
from ..preparers.standardizers import StandardizeDragos
from ..preparers.markers import PharmacophoreAtomMarker, CGRatomMarker


class Pkab(PropertyExtractor):
    def __init__(self, workpath='.', s_option=None, marker_rules=None, standardize=None, acid=True, base=True,
                 cgr_reverse=False, is_reaction=False,
                 cgr_marker=None, cgr_marker_preprocess=None, cgr_marker_postprocess=None, cgr_stereo=False):

        if is_reaction:
            if not cgr_marker:
                raise Exception('only cgr marker can work with reactions')
        elif cgr_marker:
            raise Exception('for cgr marker is_reaction should be True')

        PropertyExtractor.__init__(self, s_option)

        self.__phm_marker = PharmacophoreAtomMarker(marker_rules, workpath) if marker_rules else None

        self.__cgr_marker = CGRatomMarker(cgr_marker, preprocess=cgr_marker_preprocess,
                                          postprocess=cgr_marker_postprocess,
                                          stereo=cgr_stereo, reverse=cgr_reverse) if cgr_marker else None

        self.__dragos_std = StandardizeDragos(standardize) if standardize is not None and not is_reaction else None
        self.__workpath = workpath
        self.__reverse = cgr_reverse
        self.__acid = acid
        self.__base = base

        locs = locals()
        tmp = dict((x, locs[x]) for x in self.__optional_configs if locs[x])
        tmp.update((x, y) for x, y in (('acid', acid), ('base', base)) if not y)
        self.__config = tmp

    __optional_configs = ('s_option', 'marker_rules', 'standardize', 'cgr_reverse', 'cgr_marker',
                          'cgr_marker_preprocess', 'cgr_marker_postprocess', 'cgr_stereo', 'is_reaction')

    def get_config(self):
        return self.__config

    def set_work_path(self, workpath):
        self.__workpath = workpath
        if self.__phm_marker:
            self.__phm_marker.set_work_path(workpath)

    def get(self, structures, **_):
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

        prop, doubles, used_str = [], [], []

        p = Popen([CXCALC] + 'pka -x 50 -i -50 -a 8 -b 8 -P dynamic -m micro'.split(),
                  stdout=PIPE, stdin=PIPE, stderr=STDOUT)

        with StringIO() as f:
            writer = SDFwrite(f)
            for s_numb, s in enumerate(structures):
                if self.__cgr_marker:
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
                elif self.__phm_marker:
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
            if self.__cgr_marker or self.__phm_marker:
                old_k = None
                for lk, v in zip(doubles, pk):
                    for k in (lk if self.__phm_marker else [lk]):
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

            x = DataFrame([x[1] for x in new_doubles],
                          columns=(['pka'] if self.__acid else []) + (['pkb'] if self.__base else []))

            res = dict(X=x, AD=x.isnull().any(axis=1) ^ True, Y=Series(prop, name='Property'),
                       structures=used_str)  # todo: prepare structures

            if self.__cgr_marker or self.__phm_marker:
                i = MultiIndex.from_tuples([x[0] for x in new_doubles], names=['structure', 'c.0', 'c.1'])
            else:
                i = Index(doubles, name='structure')

            res['X'].index = res['AD'].index = res['Y'].index = i
            return res

        return False
