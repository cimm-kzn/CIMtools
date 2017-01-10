# -*- coding: utf-8 -*-
#
#  Copyright 2015, 2016 Ramil Nugmanov <stsouko@live.ru>
#  This file is part of MODtools.
#
#  MODtools is free software; you can redistribute it and/or modify
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
from os import path

CHEMAXON = "https://cimm.kpfu.ru/webservices"
JCHEMBIN = '/opt/JChem/bin'
FRAGMENTOR = '/opt/fragmentor/fragmentor'
EED = '/opt/dragos/eedstart.sh'
COLOR = '/opt/dragos/colorstart.sh'
GACONF = '/opt/dragos/gaconfstarter.sh'

config_list = ('CHEMAXON', 'JCHEMBIN', 'FRAGMENTOR', 'EED', 'COLOR', 'GACONF')

if not path.exists(path.join(path.dirname(__file__), "config.ini")):
    with open(path.join(path.dirname(__file__), "config.ini"), 'w') as f:
        f.write('\n'.join('%s = %s' % (x, y) for x, y in globals().items() if x in config_list))

with open(path.join(path.dirname(__file__), "config.ini")) as f:
    for line in f:
        try:
            k, v = line.split('=')
            k = k.strip()
            v = v.strip()
            if k in config_list:
                globals()[k] = int(v) if v.isdigit() else v
        except:
            pass


MOLCONVERT = path.join(JCHEMBIN, 'molconvert')
STANDARDIZER = path.join(JCHEMBIN, 'standardize')
CXCALC = path.join(JCHEMBIN, 'cxcalc')
REACTOR = path.join(JCHEMBIN, 'react')
JCSEARCH = path.join(JCHEMBIN, 'jcsearch')
PMAPPER = path.join(JCHEMBIN, 'pmapper')
