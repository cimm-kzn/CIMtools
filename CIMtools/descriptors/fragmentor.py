# -*- coding: utf-8 -*-
#
#  Copyright 2015-2017 Ramil Nugmanov <stsouko@live.ru>
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
from CGRtools.CGRpreparer import CGRcombo
from CGRtools.files.SDFrw import SDFwrite
from itertools import tee
from os import close, remove
from os.path import join, exists, devnull
from pandas import DataFrame, Series
from shutil import rmtree
from subprocess import call
from sys import stderr
from tempfile import mkdtemp, mkstemp
from .basegenerator import BaseGenerator
from ..config import FRAGMENTOR
from ..preparers.colorize import Colorize
from ..preparers.markers import PharmacophoreAtomMarker, CGRatomMarker
from ..preparers.standardizers import StandardizeDragos


class OpenFiles(object):
    def __init__(self, files, flags):
        if isinstance(files, str):
            files = [files]
        if isinstance(flags, str):
            flags = [flags]

        if len(flags) == 1:
            flags = flags * len(files)
        elif len(flags) != len(files):
            raise Exception('number of flags should be equal to number of files or equal to one')

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
                 s_option=None, fragment_type=3, min_length=2, max_length=10, colorname=None, marked_atom=0,
                 cgr_dynbonds=0, xml=None, doallways=False, useformalcharge=False, atompairs=False,
                 fragmentstrict=False, getatomfragment=False, overwrite=True, header=None,
                 marker_rules=None, standardize=False, docolor=False,
                 cgr_marker=None, cgr_marker_preprocess=None, cgr_marker_postprocess=None, cgr_reverse=False,
                 cgr_type=None, cgr_extralabels=False, cgr_b_templates=None, cgr_m_templates=None,
                 cgr_isotope=False, cgr_element=True, cgr_stereo=False, is_reaction=False):
        """
        Fragmentor wrapper
        :param workpath: path for temp files.
        :param version: fragmentor version
        :param s_option: modeling attribute
        :param fragment_type: 
        :param min_length: 
        :param max_length: 
        :param colorname: 
        :param marked_atom: 
        :param cgr_dynbonds: 
        :param xml: 
        :param doallways: 
        :param useformalcharge: 
        :param atompairs: 
        :param fragmentstrict: 
        :param getatomfragment: 
        :param overwrite: 
        :param header: 
        :param marker_rules: Dragos atom marker procedure. For molecules only. string with chemaxon Pmapper rules xml.
        :param standardize: Dragos standardization procedure. For molecules only.
          string with chemaxon standardizer rules (xml or ..- separated params)
          if False - skipped. if None - use predefined rules.
        :param docolor: Automatic coloring utility. see colorname.
          string with chemaxon standardizer rules (xml or ..- separated params)
          if False - skipped. if None - use predefined rules.
        :param cgr_type: type of generated CGR. see CGRtools help.
          if None - skipped. usable with is_reaction=True only.
        :param cgr_b_templates: list of reactioncontainers with termplates for autobalancing reactions
        :param cgr_m_templates: list of reactioncontainers with termplates for reactions mapping correction
        :param cgr_extralabels: match neighbors and hyb marks in substructure search. 
          Need for CGRatomMarker or CGRcombo balanser/remapper
        :param cgr_isotope: match isotope. see cgr_extralabels
        :param cgr_element: match stereo. see cgr_extralabels
        :param cgr_stereo: match stereo.
        :param cgr_reverse: for reactions only. use product structure for calculation
        :param cgr_marker: CGRtools.CGRreactor templates  
        :param cgr_marker_preprocess: prepare (if need) reaction structure for marker 
          string with chemaxon standardizer rules (xml or ..- separated params)
        :param cgr_marker_postprocess: postprocess (if need) after marker [previously prepared] reaction structure  
          string with chemaxon standardizer rules (xml or ..- separated params)
        :param is_reaction: True for reaction data
        """

        if is_reaction:
            if not (cgr_type or cgr_marker):
                raise Exception('only cgr or cgr marker can work with reactions')
            if standardize or standardize is None:
                raise Exception('standardize can work only with molecules')
            if marker_rules:
                raise Exception('pharmacophore atom marker can work only with molecules')
        elif cgr_type or cgr_marker:
            raise Exception('for cgr or cgr marker is_reaction should be True')

        self.__is_reaction = is_reaction
        self.__workpath = workpath

        BaseGenerator.__init__(self, s_option=s_option)

        self.__preprocess = any([marker_rules, standardize, standardize is None, cgr_type, cgr_marker, docolor])

        if cgr_type:
            self.__cgr = CGRcombo(cgr_type=cgr_type, extralabels=cgr_extralabels, isotope=cgr_isotope,
                                  element=cgr_element, stereo=cgr_stereo,
                                  b_templates=cgr_b_templates, m_templates=cgr_m_templates)

        if marker_rules:
            self.__marker = PharmacophoreAtomMarker(marker_rules, workpath)
        elif cgr_marker:
            self.__marker = CGRatomMarker(cgr_marker, preprocess=cgr_marker_preprocess,
                                          postprocess=cgr_marker_postprocess, extralabels=cgr_extralabels,
                                          isotope=cgr_isotope, element=cgr_element, stereo=cgr_stereo,
                                          b_templates=cgr_b_templates, m_templates=cgr_m_templates, reverse=cgr_reverse)

        self.markers = (self.__marker.get_count() if self.__marker is not None else None)
        self.__work_files = self.markers or 1

        if standardize or standardize is None:
            self.__dragos_std = StandardizeDragos(rules=standardize)

        if docolor or docolor is None:
            self.__do_color = Colorize(docolor, workpath)

        self.__frag_version = ('-%s' % version) if version else ''

        self.__head_dump = {}
        self.__head_size = {}
        self.__head_dict = {}
        self.__head_cols = {}
        self.__head_exec = {}

        if header:
            self.__gen_header = False
            self.__manual_header = True
            headers = header if isinstance(header, list) else [header]
            if len(headers) != self.__work_files:
                raise Exception('number header files should be equal to number of markers or 1')

            for n, h in enumerate(headers):
                (self.__head_dump[n], self.__head_dict[n], self.__head_cols[n],
                 self.__head_size[n]) = self.__parse_header(h)
            self.__prepare_headers()

        elif header is not None:
            self.__gen_header = False
            self.__headerless = True

        tmp = ['-f', 'SVM', '-t', str(fragment_type), '-l', str(min_length), '-u', str(max_length)]

        if colorname:
            tmp.extend(['-c', colorname])
        if marked_atom:
            tmp.extend(['-m', str(marked_atom)])
        if cgr_dynbonds:
            tmp.extend(['-d', str(cgr_dynbonds)])
        if xml:
            tmp.extend(['-x', xml])
        if doallways:
            tmp.append('--DoAllWays')
        if atompairs:
            tmp.append('--AtomPairs')
        if useformalcharge:
            tmp.append('--UseFormalCharge')
        if fragmentstrict:
            tmp.append('--StrictFrg')
        if getatomfragment:
            tmp.append('--GetAtomFragment')
        if not overwrite:
            tmp.append('--Pipe')

        self.__exec_params = tmp

        locs = locals()
        tmp = dict(fragment_type=fragment_type, min_length=min_length, max_length=max_length)
        tmp.update((x, y) for x, y in (('overwrite', overwrite), ('cgr_element', cgr_element)) if not y)
        tmp.update((x, locs[x]) for x in self.__optional_configs if locs[x])
        self.__config = tmp

    __optional_configs = ('version', 's_option', 'colorname', 'marked_atom', 'cgr_dynbonds', 'xml', 'doallways',
                          'useformalcharge',  'atompairs', 'fragmentstrict', 'getatomfragment', 'header',
                          'marker_rules', 'standardize', 'docolor', 'cgr_marker', 'cgr_marker_preprocess',
                          'cgr_marker_postprocess', 'cgr_reverse', 'cgr_type', 'cgr_extralabels', 'cgr_b_templates',
                          'cgr_m_templates', 'cgr_isotope', 'cgr_stereo', 'is_reaction')
    __gen_header = True
    __manual_header = False
    __headerless = False
    __marker = None
    __dragos_std = None
    __do_color = None
    __cgr = None

    def get_config(self):
        return self.__config

    def flush(self):
        if not (self.__gen_header or self.__manual_header or self.__headerless):
            self.__gen_header = True
            self.__head_exec = {}
            self.__head_dump = {}
            self.__head_size = {}
            self.__head_dict = {}
            self.__head_cols = {}

    @property
    def __fragmentor(self):
        return '%s%s' % (FRAGMENTOR, self.__frag_version)

    @staticmethod
    def __parse_header(header):
        with open(header, encoding='utf-8') as f:
            head_dump = f.read()
            head_dict = {int(k[:-1]): v for k, v in (i.split() for i in head_dump.splitlines())}
            head_columns = list(head_dict.values())
            head_size = len(head_dict)
        return head_dump, head_dict, head_columns, head_size

    def set_work_path(self, workpath):
        self.delete_work_path()
        if hasattr(self.__marker, 'set_work_path'):
            self.__marker.set_work_path(workpath)
        if self.__do_color:
            self.__do_color.set_work_path(workpath)

        self.__workpath = workpath
        if not (self.__headerless or self.__gen_header):
            self.__prepare_headers()

    def delete_work_path(self):
        if hasattr(self.__marker, 'delete_work_path'):
            self.__marker.delete_work_path()
        if self.__do_color:
            self.__do_color.delete_work_path()

        if not (self.__headerless or self.__gen_header) and self.__head_exec:
            for n in range(self.__work_files):
                remove(self.__head_exec.pop(n))

    def pickle(self):
        if hasattr(self.__marker, 'pickle'):
            self.__marker.pickle()

    def __prepare_headers(self):
        for n in range(self.__work_files):
            fd, header = mkstemp(prefix='frg_', suffix='.hdr', dir=self.__workpath)
            with open(header, 'w', encoding='utf-8') as f:
                f.write(self.__head_dump[n])
            close(fd)

            self.__head_exec[n] = header

    def prepare(self, structures, **_):
        """ PMAPPER and Standardizer works only with molecules. NOT CGR!
        :param structures: list of MoleculeContainers or ReactionContainers (work only in CGR or CGR-marked atoms mode)
        """
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

            elif self.__marker:
                structures = self.__marker.get(structures)

            if not structures:
                return False

        work_dir = mkdtemp(prefix='frg_', dir=self.__workpath)
        work_sdf_files = [join(work_dir, "frg_%d.sdf" % x) for x in range(self.__work_files)]
        work_files = [join(work_dir, "frg_%d" % x) for x in range(self.__work_files)]

        with OpenFiles(work_sdf_files, 'w') as f:
            writers = [SDFwrite(x) for x in f]
            prop, doubles, used_str = self.write_prepared(structures, writers)

        tx, td = [], []
        for n, (work_file, work_sdf) in enumerate(zip(work_files, work_sdf_files)):
            work_file_svm = '%s.svm' % work_file
            work_file_hdr = '%s.hdr' % work_file
            execparams = [self.__fragmentor, '-i', work_sdf, '-o', work_file]
            if not (self.__headerless or self.__gen_header):
                execparams.extend(['-h', self.__head_exec[n]])
            execparams.extend(self.__exec_params)

            print(' '.join(execparams), file=stderr)
            with open(devnull, 'w') as silent:
                exitcode = call(execparams, stdout=silent, stderr=silent) == 0

            if exitcode and exists(work_file_svm) and exists(work_file_hdr):
                if self.__headerless:
                    head_dict, head_cols, head_size = self.__parse_header(work_file_hdr)[1:]
                elif self.__gen_header:  # dump header if don't set on first run
                    self.__head_dump[n], self.__head_dict[n], self.__head_cols[n], self.__head_size[n] = \
                        _, head_dict, head_cols, head_size = self.__parse_header(work_file_hdr)

                    if n + 1 == len(work_files):  # disable header generation
                        self.__gen_header = False
                        self.__prepare_headers()
                else:
                    head_dict, head_cols, head_size = self.__head_dict[n], self.__head_cols[n], self.__head_size[n]

                x, d = self.__parse_fragmentor_output(work_file_svm, head_dict, head_cols, head_size)
                tx.append(x)
                td.append(d)
            else:
                rmtree(work_dir)
                raise Exception('Fragmentor execution FAILED')

        rmtree(work_dir)
        return tx, prop, td, doubles, used_str

    @staticmethod
    def __parse_fragmentor_output(svm_file, head_dict, head_cols, head_size):
        vector, ad = [], []
        with open(svm_file) as sf:
            for frag in sf:
                _, *x = frag.split()
                ad.append(True)
                tmp = {}  # X vector
                for i in x:
                    k, v = (int(x) for x in i.split(':'))
                    if k <= head_size:
                        tmp[head_dict[k]] = v
                    elif v != 0:
                        ad[-1] = False
                        break
                vector.append(tmp)

        return DataFrame(vector, columns=head_cols).fillna(0), Series(ad)
