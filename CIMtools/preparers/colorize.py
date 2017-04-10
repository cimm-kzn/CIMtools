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
from os import close, remove
from os.path import join, dirname
from shutil import rmtree
from subprocess import call
from tempfile import mkdtemp, mkstemp
from ..config import COLOR


class Colorize(object):
    def __init__(self, standardize, workpath='.'):
        self.__standardize = self.__dump_rules(standardize)
        self.set_work_path(workpath)

    __std_file = None

    @staticmethod
    def __dump_rules(rules):
        with rules or open(join(dirname(__file__), "standardrules_dragos.rules")) as f:
            rules = f.read()
        return rules

    def set_work_path(self, workpath):
        self.delete_work_path()

        self.__workpath = workpath
        fd, self.__std_file = mkstemp(prefix='clr_', suffix='.xml', dir=workpath)
        with open(self.__std_file, 'w') as f:
            f.write(self.__standardize)
        close(fd)

    def delete_work_path(self):
        if self.__std_file is not None:
            remove(self.__std_file)
            self.__std_file = None

    def get(self, structure):
        work_dir = mkdtemp(prefix='clr_', dir=self.__workpath)
        input_file = join(work_dir, 'colorin.sdf')
        out_file = join(work_dir, 'colorout.sdf')

        with open(input_file, 'w') as f:
            out = SDFwrite(f)
            for i in (structure if isinstance(structure, list) else [structure]):
                out.write(i)

        res = None
        if call([COLOR, input_file, out_file, self.__std_file]) == 0:
            with open(out_file) as f:
                res = SDFread(f, remap=False).read()

        rmtree(work_dir)
        return res or False
