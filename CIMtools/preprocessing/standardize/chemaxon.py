# -*- coding: utf-8 -*-
#
#  Copyright 2018 Ramil Nugmanov <stsouko@live.ru>
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
from CGRtools.files import MRVread, MRVwrite
from io import StringIO, BytesIO
from os import close, getenv
from pathlib import Path
from requests import post
from requests.exceptions import RequestException
from subprocess import run, PIPE
from sklearn.base import BaseEstimator
from tempfile import mkstemp
from ...base import CIMtoolsTransformerMixin
from ...exceptions import ConfigurationError
from ...utils import iter2array


CHEMAXON_REST = getenv('CHEMAXON_REST')


class StandardizeChemAxon(BaseEstimator, CIMtoolsTransformerMixin):
    def __init__(self, rules, workpath='.'):
        self.rules = rules
        self.set_work_path(workpath)

    def __getstate__(self):
        return {k: v for k, v in super().__getstate__().items()
                if not k.startswith('_StandardizeChemAxon__') and k != 'workpath'}

    def __setstate__(self, state):
        super().__setstate__(state)
        self.set_work_path('.')

    def __del__(self):
        self.delete_work_path()

    def set_params(self, **params):
        if params:
            super().set_params(**params)
            self.set_work_path(self.workpath)
        return self

    def set_work_path(self, workpath):
        self.workpath = workpath
        self.delete_work_path()

        fd, fn = mkstemp(prefix='std_', suffix='.xml', dir=workpath)
        self.__config = Path(fn)
        with self.__config.open('w') as f:
            f.write(self.rules)
        close(fd)

    def delete_work_path(self):
        if self.__config is not None:
            self.__config.unlink()
            self.__config = None

    def transform(self, x):
        x = super().transform(x)
        return iter2array(self.__processor_m(x) if x.size > 1 or not CHEMAXON_REST else self.__processor_s(x[0]),
                          allow_none=True)

    def __processor_m(self, structures):
        with StringIO() as f:
            with MRVwrite(f) as w:
                for s in structures:
                    w.write(s)
            tmp = f.getvalue().encode()
        try:
            p = run(['standardize', '-c', str(self.__config), '-f', 'mrv', '-g'], input=tmp, stdout=PIPE, stderr=PIPE)
        except FileNotFoundError as e:
            raise ConfigurationError from e

        if p.returncode != 0:
            raise ConfigurationError(p.stderr.decode())

        with BytesIO(p.stdout) as f, MRVread(f) as r:
            res = r.read()
            for x in p.stderr.decode().split('\n'):
                if x.startswith('SEVERE: Error at molecule No.'):
                    res.insert(int(x[29:].split(':', 1)[0]) - 1, None)
            return res

    def __processor_s(self, structure):
        with StringIO() as f:
            with MRVwrite(f) as w:
                w.write(structure)
            data = dict(structure=f.getvalue(), parameters='mrv',
                        filterChain=[dict(filter='standardizer',
                                          parameters=dict(standardizerDefinition=self.rules))])
        try:
            q = post(f'{CHEMAXON_REST}/rest-v0/util/calculate/molExport', json=data, timeout=20)
        except RequestException as e:
            raise ConfigurationError from e

        if q.status_code not in (201, 200):
            return [None]

        res = q.json()
        if not res:
            return [None]

        with BytesIO(res['structure'].encode()) as f, MRVread(f) as r:
            return r.read()

    __config = None


__all__ = ['StandardizeChemAxon']
