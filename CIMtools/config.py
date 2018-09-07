# -*- coding: utf-8 -*-
#
#  Copyright 2015-2018 Ramil Nugmanov <stsouko@live.ru>
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
from pathlib import Path
from sys import stderr
from traceback import format_exc


CHEMAXON_REST = None

config_list = 'CHEMAXON_REST',
config_dirs = [Path('~/.CIMtools.ini').expanduser(), Path('/etc/CIMtools.ini')]

if not any(x.exists() for x in config_dirs):
    with config_dirs[0].open('w') as f:
        f.write('#CHEMAXON_REST=http://url/webservices')


with next(x for x in config_dirs if x.exists()).open() as f:
    for n, line in enumerate(f, start=1):
        try:
            line = line.strip()
            if line and not line.startswith('#'):
                k, v = line.split('=')
                k = k.rstrip()
                v = v.lstrip()
                if k in config_list:
                    globals()[k] = int(v) if v.isdigit() else v == 'True' if v in ('True', 'False', '') else v
        except ValueError:
            print('line %d\n\n%s\n consist errors: %s' % (n, line, format_exc()), file=stderr)
