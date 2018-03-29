# -*- coding: utf-8 -*-
#
#  Copyright 2017, 2018 Ramil Nugmanov <stsouko@live.ru>
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
from CGRtools.containers import MoleculeContainer
from CGRtools.files import SDFread, SDFwrite
from os import close
from pathlib import Path
from shutil import rmtree
from sklearn.base import BaseEstimator
from subprocess import run, PIPE
from tempfile import mkstemp, mkdtemp
from .common import iter2array, TransformerMixin
from ..config import COLOR
from ..exceptions import ConfigurationError


class Colorize(BaseEstimator, TransformerMixin):
    def __init__(self, standardize=None, workpath='.'):
        self.standardize = standardize or self.__load_rules()
        self.set_work_path(workpath)

    def __getstate__(self):
        return {k: v for k, v in super().__getstate__().items() if not k.startswith('_Colorize__')}

    def __setstate__(self, state):
        super().__setstate__(state)
        self.set_work_path('.')

    def __del__(self):
        self.delete_work_path()

    def get_params(self, *args, **kwargs):
        return {k: v for k, v in super().get_params(*args, **kwargs).items() if k != 'workpath'}

    def set_params(self, **params):
        super().set_params(**params)
        if params:
            self.set_work_path(str(self.__config.parent))
        return self

    @staticmethod
    def __load_rules():
        with (Path(__file__).parent / 'standardize' / 'horvat.xml').open() as f:
            out = f.read()
        return out

    def set_work_path(self, workpath):
        self.delete_work_path()

        fd, fn = mkstemp(prefix='clr_', suffix='.xml', dir=workpath)
        self.__config = Path(fn)
        with self.__config.open('w') as f:
            f.write(self.standardize)
        close(fd)

    def delete_work_path(self):
        if self.__config is not None:
            self.__config.unlink()
            self.__config = None

    def transform(self, x):
        x = super().transform(x)

        work_dir = Path(mkdtemp(prefix='clr_', dir=str(self.__config.parent)))
        input_file = work_dir / 'colorin.sdf'
        out_file = work_dir / 'colorout.sdf'

        with input_file.open('w') as f, SDFwrite(f) as w:
            for s in x:
                w.write(s)

        try:
            p = run([COLOR, str(input_file), str(out_file), str(self.__config)], stderr=PIPE, stdout=PIPE)
        except FileNotFoundError as e:
            raise ConfigurationError(e)

        if p.returncode != 0:
            raise ConfigurationError(p.stderr.decode())

        with out_file.open() as f:
            res = SDFread(f, remap=False).read()

        # todo: implement error structure handler

        rmtree(str(work_dir))
        return iter2array(res, allow_none=True)

    __config = None
    _dtype = MoleculeContainer
