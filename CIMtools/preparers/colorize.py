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
from os import remove
from os.path import join, exists, dirname
from subprocess import call
from CGRtools.files.SDFrw import SDFread, SDFwrite
from ..config import COLOR


class Colorize(object):
    def __init__(self, standardize, workpath=None):
        self.__standardize = self.__dump_rules(standardize)
        self.__input_file = './colorin.sdf'
        self.__out_file = './colorout.sdf'
        self.__std_file = './colorstd.xml'
        if workpath is not None:
            self.set_work_path(workpath)
        else:
            self.__load_rules()

    @staticmethod
    def __dump_rules(rules):
        with rules or open(join(dirname(__file__), "standardrules_dragos.rules")) as f:
            rules = f.read()
        return rules

    def __load_rules(self):
        with open(self.__std_file, 'w') as f:
            f.write(self.__standardize)

    def set_work_path(self, workpath):
        self.__input_file = join(workpath, 'colorin.sdf')
        self.__out_file = join(workpath, 'colorout.sdf')
        self.__std_file = join(workpath, 'colorstd.xml')
        self.__load_rules()

    def get(self, structure):
        if exists(self.__out_file):
            remove(self.__out_file)
        with open(self.__input_file, 'w') as f:
            out = SDFwrite(f)
            for i in (structure if isinstance(structure, list) else [structure]):
                out.write(i)

        if call([COLOR, self.__input_file, self.__out_file, self.__std_file]) == 0:
            with open(self.__out_file) as f:
                res = SDFread(f, remap=False).read()
                if res:
                    return res
        return False
