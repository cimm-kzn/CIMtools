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
from CGRtools.containers import MoleculeContainer
from CGRtools.files import SDFwrite
from io import StringIO
from itertools import count
from os import close
from pathlib import Path
from sklearn.base import BaseEstimator
from subprocess import run, PIPE
from tempfile import mkstemp
from ..common import iter2array, TransformerMixin
from ...exceptions import ConfigurationError


class AtomMarkerPharmacophore(BaseEstimator, TransformerMixin):
    def __init__(self, marker_rules, workpath='.'):
        self.marker_rules = marker_rules
        self.__init()
        self.set_work_path(workpath)

    def __init(self):
        self._marks = {}
        self.__marks_counter = count(1)

    def __getstate__(self):
        return {k: v for k, v in super().__getstate__().items()
                if not k.startswith('_AtomMarkerPharmacophore__') and k != 'workpath'}

    def __setstate__(self, state):
        super().__setstate__(state)
        self.__marks_counter = count(len(self._marks) + 1)
        self.set_work_path('.')

    def __del__(self):
        self.delete_work_path()

    def set_params(self, **params):
        if params:
            super().set_params(**params)
            self.__init()
            self.set_work_path(self.workpath)
        return self

    def set_work_path(self, workpath):
        self.workpath = workpath
        self.delete_work_path()

        fd, fn = mkstemp(prefix='pmp_', suffix='.xml', dir=workpath)
        self.__config = Path(fn)
        with self.__config.open('w') as f:
            f.write(self.marker_rules)
        close(fd)

    def delete_work_path(self):
        if self.__config is not None:
            self.__config.unlink()
            self.__config = None

    def transform(self, x):
        """
        marks atoms in 7th col of sdf.
        """
        x = super().transform(x)

        with StringIO() as f, SDFwrite(f) as w:
            for s in x:
                w.write(s)
            tmp = f.getvalue().encode()

        try:
            p = run(['pmapper', '-g', '-c', str(self.__config)], input=tmp, stdout=PIPE, stderr=PIPE)
        except FileNotFoundError as e:
            raise ConfigurationError(e)

        if p.returncode != 0:
            raise ConfigurationError(p.stderr.decode())

        marks = p.stdout.decode().split()
        for x in p.stderr.decode().split('\n'):
            if x.startswith('Import error when reading molecule'):
                marks.insert(int(x[34:].split(':', 1)[0]) - 1, '')

        result = []
        for s, mark in zip(x, marks):
            mark = mark.split(';')
            if not any(mark):
                result.append(None)
            else:
                s = s.copy()
                for (_, a), m in zip(s.nodes(data=True), mark):
                    if m:
                        a['mark'] = self._marks.get(m) or self._marks.setdefault(m, str(next(self.__marks_counter)))
                result.append(s)

        return iter2array(result, allow_none=True)

    __config = None
    _dtype = MoleculeContainer


__all__ = ['AtomMarkerPharmacophore']
