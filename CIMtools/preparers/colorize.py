# -*- coding: utf-8 -*-
#
#  Copyright 2017 Ramil Nugmanov <stsouko@live.ru>
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
from CGRtools.files.SDFrw import SDFread, SDFwrite
from os import close
from pathlib import Path
from shutil import rmtree
from subprocess import call
from tempfile import mkdtemp, mkstemp
from ..config import COLOR


class Colorize(object):
    def __init__(self, standardize=None, workpath='.'):
        self.__standardize = standardize or self.__load_rules()
        self.set_work_path(workpath)

    def pickle(self):
        return dict(standardize=self.__standardize)

    @classmethod
    def unpickle(cls, config):
        if {'standardize'}.difference(config):
            raise Exception('Invalid config')
        return cls(config['standardize'])

    @staticmethod
    def __load_rules():
        with (Path(__file__).parent / 'standardrules_dragos.xml').open() as f:
            out = f.read().strip()
        return out

    def set_work_path(self, workpath):
        self.delete_work_path()
        self.__workpath = Path(workpath)
        fd, fn = mkstemp(prefix='clr_', suffix='.xml', dir=workpath)
        self.__std_file = Path(fn)
        with self.__std_file.open('w') as f:
            f.write(self.__standardize)
        close(fd)

    def delete_work_path(self):
        if self.__std_file is not None:
            self.__std_file.unlink()
            self.__std_file = None

    def get(self, structure):
        work_dir = Path(mkdtemp(prefix='clr_', dir=str(self.__workpath)))
        input_file = work_dir / 'colorin.sdf'
        out_file = work_dir / 'colorout.sdf'

        with open(input_file, 'w') as f:
            out = SDFwrite(f)
            for i in (structure if isinstance(structure, list) else [structure]):
                out.write(i)

        res = None
        if call([COLOR, str(input_file), str(out_file), self.__std_file]) == 0:
            with out_file.open() as f:
                res = SDFread(f, remap=False).read()

        rmtree(str(work_dir))
        return res or False

    __std_file = None
