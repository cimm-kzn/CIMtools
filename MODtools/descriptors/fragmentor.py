# -*- coding: utf-8 -*-
#
#  Copyright 2015-2017 Ramil Nugmanov <stsouko@live.ru>
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
from os.path import join, exists
from subprocess import call
from sys import stderr
from pandas import DataFrame, Series
from itertools import tee
from CGRtools.CGRpreparer import CGRcombo
from CGRtools.files.SDFrw import SDFwrite
from .basegenerator import BaseGenerator
from ..config import FRAGMENTOR
from ..preparers.colorize import Colorize
from ..preparers.standardizers import StandardizeDragos
from ..preparers.markers import PharmacophoreAtomMarker, CGRatomMarker


class OpenFiles(object):
    def __init__(self, files, flags):
        if isinstance(files, str):
            files = [files]
        if isinstance(flags, str):
            flags = [flags]
        assert len(flags) == len(files)
        self.files = files
        self.flags = flags

    def __enter__(self):
        self.fhs = []
        for f, fl in zip(self.files, self.flags):
            self.fhs.append(open(f, fl))
        return self.fhs

    def __exit__(self, type_, value, traceback):
        for f in self.fhs:
            f.close()


def pairwise(iterable):
    """"s -> (s0,s1), (s1,s2), (s2, s3), ..."""
    a, b = tee(iterable)
    next(b, None)
    return zip(a, b)


class Fragmentor(BaseGenerator):
    def __init__(self, workpath='.', version=None,
                 s_option=None, fragment_type=3, min_length=2, max_length=10, color_name=None, marked_atom=0,
                 cgr_dyn_bonds=0, xml=None, do_allways=False, use_formal_charge=False, atom_pairs=False,
                 fragment_strict=False, get_atom_fragment=False, overwrite=True, header=None,
                 marker_rules=None, standardize=None, do_color=None,
                 cgr_marker=None, cgr_marker_preprocess=None, cgr_marker_postprocess=None, cgr_reverse=False,
                 cgr_type=None, cgr_extralabels=False, cgr_b_templates=None, cgr_m_templates=None,
                 cgr_isotope=False, cgr_element=True, cgr_stereo=False, is_reaction=False):

        if is_reaction and not (cgr_type or cgr_marker):
            raise Exception('only cgr or cgr marker can work with reactions')

        self.__is_reaction = is_reaction

        BaseGenerator.__init__(self, workpath=workpath, s_option=s_option)

        self.__preprocess = any(x is not None for x in (marker_rules, standardize, cgr_type, cgr_marker, do_color))

        self.__phm_marker = PharmacophoreAtomMarker(marker_rules, workpath) if marker_rules else None

        self.__cgr = CGRcombo(cgr_type=cgr_type, extralabels=cgr_extralabels,
                              isotope=cgr_isotope, element=cgr_element, stereo=cgr_stereo,
                              b_templates=cgr_b_templates, m_templates=cgr_m_templates) if cgr_type else None

        self.__cgr_marker = CGRatomMarker(cgr_marker, preprocess=cgr_marker_preprocess,
                                          postprocess=cgr_marker_postprocess,
                                          stereo=cgr_stereo, reverse=cgr_reverse) if cgr_marker else None

        self.__dragos_std = StandardizeDragos(standardize) if standardize is not None and not is_reaction else None
        self.__do_color = Colorize(do_color, workpath) if do_color else None

        self.markers = self.__cgr_marker.get_count() if cgr_marker else \
            self.__phm_marker.get_count() if marker_rules else None

        self.__work_files = self.markers or 1

        self.__head_dump = {}
        self.__head_size = {}
        self.__head_dict = {}
        self.__head_columns = {}

        self.__frag_version = ('-%s' % version) if version else ''
        tmp = ['-f', 'SVM']

        self.__gen_header = True
        if header:
            self.__gen_header = False
            headers = header if isinstance(header, list) else [header]
            for n, h in enumerate(headers):
                self.__dump_header(n, h)
            tmp.extend(['-h', ''])

        tmp.extend(['-t', str(fragment_type), '-l', str(min_length), '-u', str(max_length)])

        if color_name:
            tmp.extend(['-c', color_name])
        if marked_atom:
            tmp.extend(['-m', str(marked_atom)])
        if cgr_dyn_bonds:
            tmp.extend(['-d', str(cgr_dyn_bonds)])
        if xml:
            tmp.extend(['-x', xml])
        if do_allways:
            tmp.append('--DoAllWays')
        if atom_pairs:
            tmp.append('--AtomPairs')
        if use_formal_charge:
            tmp.append('--UseFormalCharge')
        if fragment_strict:
            tmp.append('--StrictFrg')
        if get_atom_fragment:
            tmp.append('--GetAtomFragment')
        if not overwrite:
            tmp.append('--Pipe')

        self.__exec_params = tmp

    @property
    def __fragmentor(self):
        return '%s%s' % (FRAGMENTOR, self.__frag_version)

    def __dump_header(self, n, header):
        with header as f:
            self.__head_dump[n] = f.read()
            lines = self.__head_dump[n].splitlines()
            self.__head_size[n] = len(lines)
            self.__head_dict[n] = {int(k[:-1]): v for k, v in (i.split() for i in lines)}
            self.__head_columns[n] = list(self.__head_dict[n].values())

    def set_work_path(self, workpath):
        super(Fragmentor, self).set_work_path(workpath)
        if self.__phm_marker:
            self.__phm_marker.set_work_path(workpath)
        if self.__do_color:
            self.__do_color.set_work_path(workpath)

    def __prepare_header(self, n):
        header = join(self.workpath, "model.hdr")
        with open(header, 'w', encoding='utf-8') as f:
            f.write(self.__head_dump[n])
        self.__exec_params[self.__exec_params.index('-h') + 1] = header

    def prepare(self, structures, **_):
        """ PMAPPER and Standardizer works only with molecules. NOT CGR!
        :param structures: opened file or string io in sdf, mol or rdf, rxn formats
        rdf, rxn work only in CGR or reagent marked atoms mode
        """
        workfiles = [join(self.workpath, "frg_%d.sdf" % x) for x in range(self.__work_files)]
        outputfile = join(self.workpath, "frg")

        if self.__preprocess:
            if self.__dragos_std:
                structures = self.__dragos_std.get(structures)

            if self.__do_color:
                if self.__is_reaction:
                    for i in ('substrats', 'products'):
                        mols, shifts = [], [0]
                        for x in structures:
                            shifts.append(len(x[i]) + shifts[-1])
                            mols.extend(x[i])

                        colored = self.__do_color.get(mols)
                        if not colored:
                            return False

                        for (y, z), x in zip(pairwise(shifts), structures):
                            x[i] = colored[y: z]
                else:
                    structures = self.__do_color.get(structures)

            if not structures:
                return False

            if self.__cgr:
                structures = [self.__cgr.getCGR(x) for x in structures]

            elif self.__cgr_marker:
                structures = self.__cgr_marker.get(structures)

            elif self.__phm_marker:
                structures = self.__phm_marker.get(structures)

            if not structures:
                return False

        with OpenFiles(workfiles, ['w'] * self.__work_files) as f:
            writers = [SDFwrite(x) for x in f]
            prop, doubles, used_str = self.write_prepared(structures, writers)

        tx, td = [], []
        for n, workfile in enumerate(workfiles):
            if not self.__gen_header:
                """ prepare header if exist (normally true). run fragmentor.
                """
                self.__prepare_header(n)

            execparams = [self.__fragmentor, '-i', workfile, '-o', outputfile]
            execparams.extend(self.__exec_params)
            print(' '.join(execparams), file=stderr)
            exitcode = call(execparams) == 0

            if exitcode and exists(outputfile + '.svm') and exists(outputfile + '.hdr'):
                if self.__gen_header:  # dump header if don't set on first run
                    self.__dump_header(n, open(outputfile + '.hdr', encoding='utf-8'))
                    if n + 1 == len(workfiles):  # disable header generation
                        self.__gen_header = False
                        self.__exec_params.insert(self.__exec_params.index('-t'), '-h')
                        self.__exec_params.insert(self.__exec_params.index('-t'), '')

                x, d = self.__parse_fragmentor_output(n, outputfile)
                tx.append(x)
                td.append(d)
            else:
                return False

        return tx, prop, td, doubles, used_str

    def __parse_fragmentor_output(self, n, output_file):
        vector, ad = [], []
        with open(output_file + '.svm') as sf:
            for frag in sf:
                _, *x = frag.split()
                ad.append(True)
                tmp = {}  # X vector
                for i in x:
                    k, v = (int(x) for x in i.split(':'))
                    if k <= self.__head_size[n]:
                        tmp[self.__head_dict[n][k]] = v
                    elif v != 0:
                        ad[-1] = False
                        break
                vector.append(tmp)

        return DataFrame(vector, columns=self.__head_columns[n]).fillna(0), Series(ad)
