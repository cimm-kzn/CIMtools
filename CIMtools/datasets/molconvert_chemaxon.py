# -*- coding: utf-8 -*-
#
#  Copyright 2018 Ramil Nugmanov <stsouko@live.ru>
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
from CGRtools.files import MRVread
from io import StringIO, BytesIO
from pathlib import Path
from subprocess import run, PIPE
from ..config import MOLCONVERT
from ..exceptions import ConfigurationError
from ..preprocessing.common import iter2array


def molconvert_chemaxon(data):
    """
    molconvert wrapper
    :param data: buffer or string or path to file
    :return: array of molecules of reactions
    """
    if isinstance(data, Path):
        with data.open('rb') as f:
            data = f.read()
    elif isinstance(data, StringIO):
        data = data.read().encode()
    elif isinstance(data, BytesIO):
        data = data.read()
    elif hasattr(data, 'read'):  # check if data is open(filename, mode)
        data = data.read()
        if isinstance(data, str):
            data = data.encode()
    elif isinstance(data, str):
        data = data.encode()
    elif not isinstance(data, bytes):
        raise ValueError('invalid input')

    try:
        p = run([MOLCONVERT, '-g', 'mrv'], input=data, stdout=PIPE)
    except FileNotFoundError as e:
        raise ConfigurationError(e)

    if p.returncode != 0:
        raise ConfigurationError(p.stderr.decode())

    with BytesIO(p.stdout) as f, MRVread(f) as r:
        return iter2array(r)
