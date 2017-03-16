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
from os.path import join, exists, dirname, expanduser

CHEMAXON = "https://cimm.kpfu.ru/webservices"
JCHEM_DIR = join(expanduser('~'), 'ChemAxon/JChem')
FRAGMENTOR = join(expanduser('~'), 'fragmentor/fragmentor')

EED = 'eedstart.sh'
COLOR = 'colorstart.sh'
GACONF = 'dragosgfstarter.sh'

UTILS_DIR = expanduser('~')
GACONF_PATH = join(UTILS_DIR, 'GAconfig')
LIBSVM_PATH = join(GACONF_PATH, 'libsvm-3.20')

config_list = ('CHEMAXON', 'JCHEM_DIR', 'FRAGMENTOR')

config_save_list = ['UTILS_DIR', 'GACONF_PATH', 'LIBSVM_PATH']
config_save_list.extend(config_list)

config_dirs = [join(x, '.CIMtools.ini') for x in (dirname(__file__), expanduser('~'), '/etc')]

if not any(exists(x) for x in config_dirs):
    with open(config_dirs[1], 'w') as f:
        f.write('\n'.join('%s=%s' % (x, y or '') for x, y in globals().items() if x in config_save_list))

with open(next(x for x in config_dirs if exists(x))) as f:
    for line in f:
        try:
            k, v = line.split('=')
            k = k.strip()
            v = v.strip()
            if k in config_list:
                globals()[k] = int(v) if v.isdigit() else v == 'True' if v in ('True', 'False', '') else v
        except:
            pass

JCHEMBIN = join(JCHEM_DIR, 'bin')

MOLCONVERT = join(JCHEMBIN, 'molconvert')
STANDARDIZER = join(JCHEMBIN, 'standardize')
CXCALC = join(JCHEMBIN, 'cxcalc')
REACTOR = join(JCHEMBIN, 'react')
JCSEARCH = join(JCHEMBIN, 'jcsearch')
PMAPPER = join(JCHEMBIN, 'pmapper')
