# -*- coding: utf-8 -*-
#
#  Copyright 2015-2018 Ramil Nugmanov <stsouko@live.ru>
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
from CGRtools.containers import MoleculeContainer, QueryContainer
from CGRtools.files import SDFwrite
from logging import info
from os import close
from os.path import devnull
from pandas import DataFrame, Series
from pathlib import Path
from shutil import rmtree
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.exceptions import NotFittedError
from subprocess import call
from tempfile import mkdtemp, mkstemp
from warnings import warn
from .common import iter2array
from ..config import FRAGMENTOR
from ..exceptions import ConfigurationError


class Fragmentor(BaseEstimator, TransformerMixin):
    def __init__(self, fragment_type=3, min_length=2, max_length=10, colorname=None, marked_atom=0, cgr_dynbonds=0,
                 xml=None, doallways=False, useformalcharge=False, atompairs=False, header=None, fragmentstrict=False,
                 getatomfragment=False, overwrite=True, workpath='.', version=None, verbose=False):
        """
        ISIDA Fragmentor wrapper

        :param workpath: path for temp files.
        :param version: fragmentor version
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
        """
        self.fragment_type = fragment_type
        self.min_length = min_length
        self.max_length = max_length
        self.colorname = colorname
        self.marked_atom = marked_atom
        self.cgr_dynbonds = cgr_dynbonds
        self.xml = xml
        self.doallways = doallways
        self.useformalcharge = useformalcharge
        self.atompairs = atompairs
        self.header = header
        self.fragmentstrict = fragmentstrict
        self.getatomfragment = getatomfragment
        self.overwrite = overwrite
        self.version = version
        self.verbose = verbose

        self.__init_header(header)
        self.set_work_path(workpath)

    def __getstate__(self):
        return {k: v for k, v in super().__getstate__().items()
                if k != 'header' and
                (not k.startswith('_Fragmentor__') or
                 k in ('_Fragmentor__head_dump', '_Fragmentor__head_less', '_Fragmentor__head_generate'))}

    def __setstate__(self, state):
        super().__setstate__({k: v for k, v in state.items() if k != '_Fragmentor__head_dump'})
        if state.get('_Fragmentor__head_dump'):
            self.__load_header(state['_Fragmentor__head_dump'])
        self.set_work_path('.')

    def __del__(self):
        self.delete_work_path()

    def get_params(self, *args, **kwargs):
        init = {k: v for k, v in super().get_params(*args, **kwargs).items() if k not in ('workpath', 'header')}
        init.update(__head_generate=self.__head_generate, __head_less=self.__head_less, __head_dump=self.__head_dump)
        return init

    def set_params(self, **params):
        if params:
            super().set_params(**{k: v for k, v in params.items() if not k.startswith('__head_')})
            if '__head_dump' in params:
                try:
                    dump = params['__head_dump']
                    self.__head_generate = params['__head_generate']
                    self.__head_less = params['__head_less']
                except KeyError as e:
                    raise ConfigurationError(e)

                if dump:
                    self.__load_header(dump)
            else:
                self.__init_header(params.get('header'))

            self.set_work_path(params.get('workpath') or self.__workpath)
        return self

    def set_work_path(self, workpath):
        self.__workpath = Path(workpath)
        if self.__head_dump is not None:
            self.__prepare_header()

    def delete_work_path(self):
        if self.__head_exec is not None:
            self.__head_exec.unlink()
            self.__head_exec = None

    def finalize(self):
        """
        finalize partial fitting procedure
        """
        if not self.__head_exec:
            raise NotFittedError('fragmentor instance is not fitted yet')
        if self.__head_generate:
            self.__head_generate = False

    def _reset(self):
        """Reset internal data-dependent state.
        __init__ parameters are not touched.
        """
        if not self.__head_less:
            if not self.__head_generate:
                self.__head_generate = True
            if self.__head_dump is not None:
                self.__head_dump = self.__head_dict = self.__head_cols = self.__head_size = None

            self.delete_work_path()

    def fit(self, x, y=None):
        """Compute the header.
        """
        x = iter2array(x, dtype=(MoleculeContainer, QueryContainer))

        if self.__head_less:
            warn('Fragmentor configured to head less mode. fit unusable')
            return self

        self._reset()
        self.__prepare(x)
        return self

    def partial_fit(self, x, y=None):
        x = iter2array(x, dtype=(MoleculeContainer, QueryContainer))

        if self.__head_less:
            warn('Fragmentor configured to head less mode. fit unusable')
            return self
        if not self.__head_generate:
            raise AttributeError('partial fit impossible. transformer already fitted')

        self.__prepare(x, True)
        return self

    def transform(self, x, return_domain=False):
        if not (self.__head_less or self.__head_exec):
            raise NotFittedError('fragmentor instance is not fitted yet')

        x = iter2array(x, dtype=(MoleculeContainer, QueryContainer))
        x, d = self.__prepare(x)
        if return_domain:
            return x, d
        return x

    def fit_transform(self, x, y=None, return_domain=False):
        x = iter2array(x, dtype=(MoleculeContainer, QueryContainer))

        self._reset()
        x, d = self.__prepare(x)
        if return_domain:
            return x, d
        return x

    def __prepare(self, x, partial=False):
        work_dir = Path(mkdtemp(prefix='frg_', dir=str(self.__workpath)))
        inp_file = work_dir / 'input.sdf'
        out_file = work_dir / 'output'
        out_file_svm = work_dir / 'output.svm'
        out_file_hdr = work_dir / 'output.hdr'

        with inp_file.open('w', encoding='utf-8') as f, SDFwrite(f) as w:
            for s in x:
                w.write(s)

        execparams = self.__exec_params(inp_file, out_file)
        info(' '.join(execparams))
        if self.verbose:
            exitcode = call(execparams) == 0
        else:
            with open(devnull, 'w') as silent:
                exitcode = call(execparams, stdout=silent, stderr=silent) == 0

        if exitcode and out_file_svm.exists() and out_file_hdr.exists():
            if self.__head_less:
                head_dict, head_cols, head_size = self.__parse_header(out_file_hdr)[1:]
            else:
                if self.__head_generate:  # dump header if don't set on first run
                    self.__load_header(out_file_hdr)
                    if not partial:
                        self.__head_generate = False
                    self.__prepare_header()

                head_dict, head_cols, head_size = self.__head_dict, self.__head_cols, self.__head_size

            try:
                x, d = self.__parse_svm(out_file_svm, head_dict, head_cols, head_size)
            except Exception as e:
                raise ConfigurationError(e)
        else:
            raise ConfigurationError('Fragmentor execution FAILED')

        rmtree(str(work_dir))
        return x, d

    @staticmethod
    def __parse_svm(svm_file, head_dict, head_cols, head_size):
        vector, ad = [], []
        with svm_file.open() as sf:
            for frag in sf:
                _, *x = frag.split()
                ad.append(True)
                tmp = {}  # X vector
                for i in x:
                    k, v = i.split(':')
                    k, v = int(k), int(v)
                    if k <= head_size:
                        tmp[head_dict[k]] = v
                    elif v != 0:
                        ad[-1] = False
                        break
                vector.append(Series(tmp))

        return DataFrame(vector, columns=head_cols).fillna(0), Series(ad)

    def __exec_params(self, inp, out):
        tmp = ['%s-%s' % (FRAGMENTOR, self.version) if self.version else FRAGMENTOR, '-i', str(inp), '-o', str(out)]

        if self.__head_exec:
            tmp.extend(('-h', str(self.__head_exec)))

        tmp.extend(('-f', 'SVM', '-t', str(self.fragment_type), '-l', str(self.min_length), '-u', str(self.max_length)))

        if self.colorname:
            tmp.extend(['-c', self.colorname])
        if self.marked_atom:
            tmp.extend(['-m', str(self.marked_atom)])
        if self.cgr_dynbonds:
            tmp.extend(['-d', str(self.cgr_dynbonds)])
        if self.xml:
            tmp.extend(['-x', self.xml])
        if self.doallways:
            tmp.append('--DoAllWays')
        if self.atompairs:
            tmp.append('--AtomPairs')
        if self.useformalcharge:
            tmp.append('--UseFormalCharge')
        if self.fragmentstrict:
            tmp.append('--StrictFrg')
        if self.getatomfragment:
            tmp.append('--GetAtomFragment')
        if not self.overwrite:
            tmp.append('--Pipe')

        return tuple(tmp)

    @staticmethod
    def __parse_header(header):
        if isinstance(header, Path):
            with header.open(encoding='utf-8') as f:
                head_dump = f.read()
        else:
            head_dump = header

        head_dict = {int(k[:-1]): v for k, v in (i.split() for i in head_dump.splitlines())}
        return head_dump, head_dict, list(head_dict.values()), len(head_dict)

    def __load_header(self, header):
        self.__head_dump, self.__head_dict, self.__head_cols, self.__head_size = self.__parse_header(header)

    def __init_header(self, header):
        if header:
            if self.__head_generate:
                self.__head_generate = False
            if self.__head_less:
                self.__head_less = False
            self.__load_header(Path(header))
        elif header is not None:
            if self.__head_generate:
                self.__head_generate = False
            if not self.__head_less:
                self.__head_less = True
        else:
            if not self.__head_generate:
                self.__head_generate = True
            if self.__head_less:
                self.__head_less = False

    def __prepare_header(self):
        fd, fn = mkstemp(prefix='frg_', suffix='.hdr', dir=str(self.__workpath))

        if self.__head_exec is not None:  # remove old header
            self.__head_exec.unlink()

        self.__head_exec = header = Path(fn)
        with header.open('w', encoding='utf-8') as f:
            f.write(self.__head_dump)
        close(fd)

    __head_generate = True
    __head_less = False
    __head_dump = __head_size = __head_dict = __head_cols = __head_exec = __workpath = None
