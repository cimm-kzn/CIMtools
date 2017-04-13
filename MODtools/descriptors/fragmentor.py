# -*- coding: utf-8 -*-
#
#  Copyright 2015, 2016 Ramil Nugmanov <stsouko@live.ru>
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
import os
import subprocess as sp
import sys
import pandas as pd
from itertools import tee
from CGRtools.CGRpreparer import CGRcombo
from CGRtools.files.SDFrw import SDFwrite
from .basegenerator import BaseGenerator
from ..config import FRAGMENTOR
from ..structprepare import Pharmacophoreatommarker, StandardizeDragos, CGRatommarker, Colorize


class openFiles(object):
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
    """s -> (s0,s1), (s1,s2), (s2, s3), ..."""
    a, b = tee(iterable)
    next(b, None)
    return zip(a, b)


class Fragmentor(BaseGenerator):
    def __init__(self, workpath='.', version=None, s_option=None, fragment_type=3, min_length=2, max_length=10,
                 colorname=None, marked_atom=0, cgr_dynbonds=0, xml=None, doallways=False, useformalcharge=False,
                 atompairs=False, fragmentstrict=False, getatomfragment=False, overwrite=True, header=None,
                 marker_rules=None, standardize=None, docolor=None,
                 cgr_marker=None, cgr_marker_prepare=None, cgr_marker_postprocess=None, cgr_reverse=False,
                 cgr_type=None, cgr_extralabels=False, cgr_b_templates=None, cgr_m_templates=None,
                 cgr_isotope=False, cgr_element=True, cgr_stereo=False, is_reaction=False):

        if is_reaction and not (cgr_type or cgr_marker):
            raise Exception('only cgr or cgr marker can work with reactions')

        self.__is_reaction = is_reaction

        BaseGenerator.__init__(self, workpath=workpath, s_option=s_option)

        self.__preprocess = any(x is not None for x in (marker_rules, standardize, cgr_type, cgr_marker, docolor))

        self.__phm_marker = Pharmacophoreatommarker(marker_rules, workpath) if marker_rules else None

        self.__cgr = CGRcombo(cgr_type=cgr_type, extralabels=cgr_extralabels,
                              isotope=cgr_isotope, element=cgr_element, stereo=cgr_stereo,
                              b_templates=cgr_b_templates, m_templates=cgr_m_templates) if cgr_type else None

        self.__cgr_marker = CGRatommarker(cgr_marker, prepare=cgr_marker_prepare,
                                          postprocess=cgr_marker_postprocess,
                                          stereo=cgr_stereo, reverse=cgr_reverse) if cgr_marker else None

        self.__dragos_std = StandardizeDragos(standardize) if standardize is not None and not is_reaction else None
        self.__do_color = Colorize(docolor, workpath) if docolor else None

        self.markers = self.__cgr_marker.getcount() if cgr_marker else \
            self.__phm_marker.getcount() if marker_rules else None

        self.__workfiles = self.markers or 1

        self.__headdump = {}
        self.__headsize = {}
        self.__headdict = {}
        self.__headcolumns = {}

        self.__fragversion = ('-%s' % version) if version else ''
        tmp = ['-f', 'SVM']

        self.__genheader = True
        if header:
            self.__genheader = False
            headers = header if isinstance(header, list) else [header]
            for n, h in enumerate(headers):
                self.__dumpheader(n, h)
            tmp.extend(['-h', ''])

        tmp.extend(['-t', str(fragment_type), '-l', str(min_length), '-u', str(max_length)])

        if colorname: tmp.extend(['-c', colorname])
        if marked_atom: tmp.extend(['-m', str(marked_atom)])
        if cgr_dynbonds: tmp.extend(['-d', str(cgr_dynbonds)])
        if xml: tmp.extend(['-x', xml])
        if doallways: tmp.append('--DoAllWays')
        if atompairs: tmp.append('--AtomPairs')
        if useformalcharge: tmp.append('--UseFormalCharge')
        if fragmentstrict: tmp.append('--StrictFrg')
        if getatomfragment: tmp.append('--GetAtomFragment')
        if not overwrite: tmp.append('--Pipe')

        self.__execparams = tmp

    def __fragmentor(self):
        return '%s%s' % (FRAGMENTOR, self.__fragversion)

    def __dumpheader(self, n, header):
        with header as f:
            self.__headdump[n] = f.read()
            lines = self.__headdump[n].splitlines()
            self.__headsize[n] = len(lines)
            self.__headdict[n] = {int(k[:-1]): v for k, v in (i.split() for i in lines)}
            self.__headcolumns[n] = list(self.__headdict[n].values())

    def setworkpath(self, workpath):
        BaseGenerator.setworkpath(self, workpath)
        if self.__phm_marker:
            self.__phm_marker.setworkpath(workpath)
        if self.__do_color:
            self.__do_color.setworkpath(workpath)

    def __prepareheader(self, n):
        header = os.path.join(self.workpath, "model.hdr")
        with open(header, 'w', encoding='utf-8') as f:
            f.write(self.__headdump[n])
        self.__execparams[self.__execparams.index('-h') + 1] = header

    def prepare(self, structures, **kwargs):
        """ PMAPPER and Standardizer works only with molecules. NOT CGR!
        :param structures: opened file or string io in sdf, mol or rdf, rxn formats
        rdf, rxn work only in CGR or reagent marked atoms mode
        """
        workfiles = [os.path.join(self.workpath, "frg_%d.sdf" % x) for x in range(self.__workfiles)]
        outputfile = os.path.join(self.workpath, "frg")

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

        with openFiles(workfiles, ['w'] * self.__workfiles) as f:
            writers = [SDFwrite(x) for x in f]
            prop, doubles, used_str = self.write_prepared(structures, writers)

        tX, tD = [], []
        for n, workfile in enumerate(workfiles):
            if not self.__genheader:
                """ prepare header if exist (normally true). run fragmentor.
                """
                self.__prepareheader(n)

            execparams = [self.__fragmentor(), '-i', workfile, '-o', outputfile]
            execparams.extend(self.__execparams)
            print(' '.join(execparams), file=sys.stderr)
            exitcode = sp.call(execparams) == 0

            if exitcode and os.path.exists(outputfile + '.svm') and os.path.exists(outputfile + '.hdr'):
                if self.__genheader:  # dump header if don't set on first run
                    self.__dumpheader(n, open(outputfile + '.hdr', encoding='utf-8'))
                    if n + 1 == len(workfiles):  # disable header generation
                        self.__genheader = False
                        self.__execparams.insert(self.__execparams.index('-t'), '-h')
                        self.__execparams.insert(self.__execparams.index('-t'), '')

                X, D = self.__parsefragmentoroutput(n, outputfile)
                tX.append(X)
                tD.append(D)
            else:
                return False

        return tX, prop, tD, doubles, used_str

    def __parsefragmentoroutput(self, n, outputfile):
        vector, ad = [], []
        with open(outputfile + '.svm') as sf:
            for frag in sf:
                _, *x = frag.split()
                ad.append(True)
                tmp = {}  # X vector
                for i in x:
                    k, v = (int(x) for x in i.split(':'))
                    if k <= self.__headsize[n]:
                        tmp[self.__headdict[n][k]] = v
                    elif v != 0:
                        ad[-1] = False
                        break
                vector.append(tmp)

        return pd.DataFrame(vector, columns=self.__headcolumns[n]).fillna(0), pd.Series(ad)
