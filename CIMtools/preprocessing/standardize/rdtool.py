# -*- coding: utf-8 -*-
#
#  Copyright 2020 Zarina Ibragimova <zarinaIbr12@yandex.ru>
#  Copyright 2020 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from CGRtools.containers import ReactionContainer
from CGRtools import RDFRead, RDFWrite
from os import name
from os.path import devnull
from pandas import DataFrame
from pathlib import Path
from shutil import rmtree
from subprocess import call
from sys import prefix, exec_prefix
from tempfile import mkdtemp
from warnings import warn
from ...base import CIMtoolsTransformerMixin
from ...exceptions import ConfigurationError


class RDTool(CIMtoolsTransformerMixin):
    def __init__(self, algorithm='max', verbose=False):
        """
        :param algorithm: 'max','min','mixture'
        """
        self.algorithm = algorithm
        self.verbose = verbose

    def transform(self, x):
        x = super().transform(x)

        algorithms = ['max', 'min', 'mixture']
        if self.algorithm not in algorithms:
            raise ValueError("Invalid value for algorithm of mapping. Allowed string values are 'max','min','mixture'")
        algorithms.remove(self.algorithm)

        work_dir = Path(mkdtemp(prefix='rdt_'))
        input_file = work_dir / 're_map.rdf'
        out_folder = work_dir / 'results'

        with RDFWrite(input_file) as f:
            for num, r in enumerate(x):
                meta = r.meta.copy()
                r.meta.clear()
                r.meta['Id'] = num
                f.write(r)
                r.meta.clear()
                r.meta.update(meta)

        execparams = ['java', '-jar', jar_path, '-j', 'MAPPING', '-i', input_file, '-o', out_folder, '-rdf_id', 'Id',
                      '-' + algorithms[0], '-' + algorithms[1]]
        try:
            if self.verbose:
                exitcode = call(execparams) != 0
            else:
                with open(devnull, 'w') as silent:
                    exitcode = call(execparams, stdout=silent, stderr=silent) != 0
        except FileNotFoundError as e:
            rmtree(work_dir)
            raise ConfigurationError(e)

        if exitcode:
            rmtree(work_dir)
            raise ConfigurationError('execution failed')

        out_file = self.algorithm.upper() + '_reactions.rdf'
        x_out = RDFRead(Path(out_folder / out_file)).read()
        if len(x) != len(x_out):
            rmtree(work_dir)
            raise ValueError('invalid data')

        rmtree(work_dir)
        return DataFrame([[x] for x in x_out], columns=['reaction'])

    _dtype = ReactionContainer


def getsitepackages():
    """returns a list containing all global site-packages directories. stolen and modified site.py function
    """
    sitepackages = []
    for pr in {prefix, exec_prefix}:
        pr = Path(pr)
        if name == 'posix':
            sitepackages.append(pr / 'local/lib')
        else:
            sitepackages.append(pr)
        sitepackages.append(pr / 'lib')
    return sitepackages


for site in getsitepackages():
    jar_path = site / 'rdtool.jar'
    if jar_path.exists():
        __all__ = ['RDTool']
        break
else:
    warn('broken package installation. rdtool.jar not found', ImportWarning)
    __all__ = []
    del RDTool
