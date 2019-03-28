# -*- coding: utf-8 -*-
#
#  Copyright 2015-2019 Ramil Nugmanov <stsouko@live.ru>
#  This file is part of CIMtools.
#
#  CIMtools is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, see <https://www.gnu.org/licenses/>.
#
from CGRtools.containers import MoleculeContainer, CGRContainer
from CGRtools.files import SDFwrite
from logging import info
from os import close
from os.path import devnull
from pandas import DataFrame, Series, concat
from pathlib import Path
from shutil import rmtree
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.exceptions import NotFittedError
from subprocess import call
from tempfile import mkdtemp, mkstemp
from warnings import warn
from ..exceptions import ConfigurationError
from ..utils import iter2array


class Fragmentor(BaseEstimator, TransformerMixin):
    def __init__(self, fragment_type=3, min_length=2, max_length=10, cgr_dynbonds=0, doallways=False,
                 useformalcharge=False, header=None, workpath='.', version=None, verbose=False, remove_rare_ratio=0,
                 return_domain=False):
        """
        ISIDA Fragmentor wrapper

        :param workpath: path for temp files.
        :param version: fragmentor version. need for selecting Fragmentor executables named as fragmentor-{version}
        :param header: if None descriptors will be generated on train set
                       if False Fragmentor will work in headless mode. in this mod fit unusable and Fragmentor return
                           all found descriptors
                       else path string to existing header file acceptable
        :param remove_rare_ratio: if descriptors found on train less then given ratio it will be removed from header.
                                  if partial fit used, be sure to use finalize method.
                                  unusable if headless mode set
        :param return_domain: add AD bool column. if False molecule has new features
        """
        self.fragment_type = fragment_type
        self.min_length = min_length
        self.max_length = max_length
        self.cgr_dynbonds = cgr_dynbonds
        self.doallways = doallways
        self.useformalcharge = useformalcharge
        self.version = version
        self.verbose = verbose
        self.header = header
        self.remove_rare_ratio = remove_rare_ratio
        self.return_domain = return_domain

        self.__init_header()
        self.set_work_path(workpath)

    def __getstate__(self):
        return {k: v for k, v in super().__getstate__().items() if
                k in ('_Fragmentor__head_dump', '_Fragmentor__head_less', '_Fragmentor__head_generate') or
                k == '_Fragmentor__head_rare' and self.__head_generate or
                k not in ('header', 'workpath') and not k.startswith('_Fragmentor__')}

    def __setstate__(self, state):
        super().__setstate__({k: v for k, v in state.items() if k != '_Fragmentor__head_dump'})
        # backward compatibility with 1.4.0 - 1.4.6
        if '_Fragmentor__head_less' not in state:
            self.__head_less = False
        if '_Fragmentor__head_generate' not in state:
            self.__head_generate = True
        if 'return_domain' not in state:
            self.return_domain = False

        if state.get('_Fragmentor__head_dump'):
            self.__load_header(state['_Fragmentor__head_dump'])
        self.set_work_path('.')

    def __del__(self):
        self.delete_work_path()

    def set_params(self, **params):
        if not params:
            return self

        self._reset()
        super().set_params(**params)
        if 'header' in params:
            self.__init_header()

        self.set_work_path(self.workpath)
        return self

    def set_work_path(self, workpath):
        self.delete_work_path()
        self.workpath = workpath
        self.__workpath = Path(workpath)
        if self.__head_dict:
            self.__prepare_header()

    def delete_work_path(self):
        if self.__head_exec:
            self.__head_exec.unlink()
            self.__head_exec = None

    def finalize(self):
        """
        finalize partial fitting procedure
        """
        if self.__head_less:
            warn(f'{self.__class__.__name__} configured to head less mode. finalize unusable')
        elif not self.__head_generate:
            warn(f'{self.__class__.__name__} already finalized or fitted')
        elif not self.__head_dict:
            raise NotFittedError(f'{self.__class__.__name__} instance is not fitted yet')
        else:
            if self.remove_rare_ratio:
                self.__clean_head(*self.__head_rare)
                self.__prepare_header()
                self.__head_rare = None
            self.__head_generate = False

    def _reset(self):
        """Reset internal data-dependent state.
        __init__ parameters are not touched.
        """
        if not self.__head_less:
            if not self.__head_generate:
                self.__head_generate = True
            if self.__head_dict:
                self.__head_dump = self.__head_dict = None
            if self.__head_rare is not None:
                self.__head_rare = None

            self.delete_work_path()

    def get_feature_names(self):
        """Get feature names.

        Returns
        -------
        feature_names : list of strings
            Names of the features produced by transform.
        """
        if self.__head_less:
            raise AttributeError(f'{self.__class__.__name__} instance configured to head less mode')
        elif not self.__head_dict:
            raise NotFittedError(f'{self.__class__.__name__} instance is not fitted yet')
        return list(self.__head_dict.values())

    def fit(self, x, y=None):
        """Compute the header.
        """
        x = iter2array(x, dtype=(MoleculeContainer, CGRContainer))

        if self.__head_less:
            warn(f'{self.__class__.__name__} configured to head less mode. fit unusable')
            return self

        self._reset()
        self.__prepare(x)
        return self

    def partial_fit(self, x, y=None):
        x = iter2array(x, dtype=(MoleculeContainer, CGRContainer))

        if self.__head_less:
            warn(f'{self.__class__.__name__} configured to head less mode. fit unusable')
            return self
        if not self.__head_generate:
            raise AttributeError(f'partial fit impossible. {self.__class__.__name__} already finalized or fitted')

        self.__prepare(x, partial=True)
        return self

    def transform(self, x):
        if not (self.__head_less or self.__head_dict):
            raise NotFittedError(f'{self.__class__.__name__} instance is not fitted yet')

        x = iter2array(x, dtype=(MoleculeContainer, CGRContainer))
        x, d = self.__prepare(x, fit=False)
        if self.return_domain:
            x['AD'] = d
        return x

    def fit_transform(self, x, y=None):
        x = iter2array(x, dtype=(MoleculeContainer, CGRContainer))
        if self.__head_less:
            warn(f'{self.__class__.__name__} configured to head less mode')

        self._reset()
        x, d = self.__prepare(x, transform=True)
        if self.return_domain:
            x['AD'] = d
        return x

    @property
    def _number_of_fragments(self):
        return len(self.__head_dict)

    @property
    def _fragments(self):
        return list(self.__head_dict)

    def __prepare(self, x, partial=False, fit=True, transform=False):
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

        if not (exitcode and out_file_svm.exists() and out_file_hdr.exists()):
            raise ConfigurationError(f'{self.__class__.__name__} execution FAILED')

        if self.__head_less:
            head_dict = self.__parse_header(out_file_hdr)
        else:
            if fit:  # dump header
                self.__load_header(out_file_hdr)
                if not partial:
                    self.__head_generate = False
            head_dict = self.__head_dict

        try:
            x, d = self.__parse_svm(out_file_svm, head_dict)
        except Exception as e:
            raise ConfigurationError(e)

        rmtree(str(work_dir))

        if not self.__head_less and fit:
            if self.remove_rare_ratio:
                amount = x.astype(bool).sum()
                if partial:
                    if self.__head_rare is None:
                        self.__head_rare = (amount, len(x))
                    else:
                        self.__head_rare = (concat([amount, self.__head_rare[0]], axis=1).sum(axis=1),
                                            self.__head_rare[1] + len(x))
                else:
                    self.__clean_head(amount, len(x))

                    if transform:
                        x = x[list(self.__head_dict.values())]

            self.__prepare_header()
        return x, d

    def __clean_head(self, fragments, total):
        c = 0
        head_dict = {}
        part = (fragments / total) >= self.remove_rare_ratio
        for k, v in part.items():
            if v:
                c += 1
                head_dict[c] = k

        info('cleaned %d rare fragments' % (len(self.__head_dict) - c))
        self.__head_dict = head_dict
        self.__head_dump = self.__format_header(self.__head_dict)

    @staticmethod
    def __parse_svm(svm_file, head_dict):
        head_size = len(head_dict)
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
                vector.append(tmp)

        return DataFrame(vector, columns=list(head_dict.values())).fillna(0), Series(ad)

    def __exec_params(self, inp, out):
        tmp = [f'fragmentor-{self.version}' if self.version else 'fragmentor', '-i', str(inp), '-o', str(out)]

        if self.__head_exec:
            tmp.extend(('-h', str(self.__head_exec)))

        tmp.extend(('-f', 'SVM', '-t', str(self.fragment_type), '-l', str(self.min_length), '-u', str(self.max_length)))

        if self.cgr_dynbonds:
            tmp.extend(['-d', str(self.cgr_dynbonds)])
        if self.doallways:
            tmp.append('--DoAllWays')
        if self.useformalcharge:
            tmp.append('--UseFormalCharge')

        return tuple(tmp)

    @staticmethod
    def __format_header(head_dict):
        return '\n'.join('%d. %s' % x for x in head_dict.items())

    @staticmethod
    def __parse_header(header):
        if isinstance(header, Path):
            with header.open(encoding='utf-8') as f:
                head_dump = f.read()
        else:
            head_dump = header
        try:
            head_dict = {int(k[:-1]): v for k, v in (i.split() for i in head_dump.splitlines())}
        except ValueError as e:
            raise ConfigurationError from e
        if not head_dict:
            raise ConfigurationError('empty header')

        return head_dict

    def __load_header(self, header):
        self.__head_dict = self.__parse_header(header)
        self.__head_dump = self.__format_header(self.__head_dict)

    def __prepare_header(self):
        if not self.__head_exec:
            fd, fn = mkstemp(prefix='frg_', suffix='.hdr', dir=str(self.__workpath))
            close(fd)
            self.__head_exec = Path(fn)

        with self.__head_exec.open('w', encoding='utf-8') as f:
            f.write(self.__head_dump)

    def __init_header(self):
        header = self.header
        if header:
            self.__head_generate = False
            self.__head_less = False
            self.__load_header(Path(header))
        elif header is not None:
            self.__head_generate = False
            self.__head_less = True
        else:
            self.__head_generate = True
            self.__head_less = False

    __head_dump = __head_dict = __head_exec = __head_rare = __workpath = None


__all__ = ['Fragmentor']
