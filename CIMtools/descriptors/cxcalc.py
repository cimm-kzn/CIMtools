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
        :param acid: calc acidity
        :param base: calc basicity
        """
        if is_reaction and not cgr_marker:
            raise Exception('need CGR marker for work with reactions')

        BaseGenerator.__init__(self, workpath=workpath, s_option=s_option, marker_rules=marker_rules,
                               standardize=standardize, cgr_marker=cgr_marker, cgr_isotope=cgr_isotope,
                               cgr_marker_preprocess=cgr_marker_preprocess, cgr_extralabels=cgr_extralabels,
                               cgr_marker_postprocess=cgr_marker_postprocess, cgr_element=cgr_element,
                               cgr_stereo=cgr_stereo, cgr_b_templates=cgr_b_templates, cgr_m_templates=cgr_m_templates,
                               cgr_reverse=cgr_reverse, is_reaction=is_reaction)

        self.__init_common(cgr_reverse, acid, base)

    def _init_unpickle(self, s_option, cgr_isotope, cgr_element, cgr_stereo, cgr_marker_preprocess,
                       cgr_marker_postprocess, cgr_extralabels, cgr_reverse, is_reaction, marker_rules, cgr_marker,
                       cgr_b_templates, cgr_m_templates, standardize, acid, base, **config):
        BaseGenerator._init_unpickle(self, s_option, cgr_isotope, cgr_element, cgr_stereo, cgr_marker_preprocess,
                                     cgr_marker_postprocess, cgr_extralabels, cgr_reverse, is_reaction,
                                     marker_rules, cgr_marker, cgr_b_templates, cgr_m_templates, standardize, **config)

        self.__init_common(cgr_reverse, acid, base)

    def __init_common(self, cgr_reverse, acid, base):
        self.__reverse = cgr_reverse
        self.__acid = acid
        self.__base = base

        tmp = []
        if acid:
            tmp.append('pka')
        if base:
            tmp.append('pkb')
        self.__columns = tmp

    def pickle(self):
        return dict(acid=self.__acid, base=self.__base, **super(Pkab, self).pickle())

    @classmethod
    def unpickle(cls, config):
        args = {'acid', 'base'}
        if args.difference(config):
            raise Exception('Invalid config')
        BaseGenerator.unpickle(config)
        obj = cls.__new__(cls)
        obj._init_unpickle(**config)
        return obj

    def set_work_path(self, workpath):
        super(Pkab, self).set_work_path()

    def delete_work_path(self):
        super(Pkab, self).delete_work_path()

    def _prepare(self, structures, **_):
        if self._dragos_std:
            structures = self._dragos_std.get(structures)

        if not structures:
            return False

        if self._marker:
            structures = self._marker.get(structures)

        if not structures:
            return False

        prop, doubles = [], []

        p = Popen([CXCALC] + 'pka -x 50 -i -50 -a 8 -b 8 -P dynamic -m micro'.split(),
                  stdout=PIPE, stdin=PIPE, stderr=STDOUT)

        with StringIO() as f:
            writer = SDFwrite(f)
            for s_numb, s in enumerate(structures):
                if isinstance(self._marker, CGRatomMarker):
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
                elif isinstance(self._marker, PharmacophoreAtomMarker):
                    writer.write(s[0][0][1])
                    prop.extend([self.get_property(d[0][1].meta, marks=[x[0] for x in d]) for d in s])
                    doubles.append([([s_numb] + [x[0] for x in d]) for d in s])
                else:
                    writer.write(s)
                    prop.append(self.get_property(s.meta))
                    doubles.append(s_numb)

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
        if self._marker:
            old_k = None
            for lk, v in zip(doubles, pk):
                for k in (lk if isinstance(self._marker, PharmacophoreAtomMarker) else [lk]):
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
        i = [x[0] for x in new_doubles] if self._markers_count else doubles

        return [x], prop, [ad], i

    def _write_prepared(self, *_):
        raise Exception('Disabled method')
