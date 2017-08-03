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
from pandas.compat import BytesIO, cPickle as pkl, pickle_compat as pc, PY3


def to_pickle(obj):
    return pkl.dumps(obj, protocol=pkl.HIGHEST_PROTOCOL)


def read_pickle(dump):
    """
    Load pickled pandas object (or any other pickled object) from the specified
    bytestring

    Warning: Loading pickled data received from untrusted sources can be
    unsafe. See: http://docs.python.org/2.7/library/pickle.html

    Parameters
    ----------
    dump : bytes

    Returns
    -------
    unpickled : type of object stored in file
    """

    def try_read(data, encoding=None):
        # try with cPickle
        # try with current pickle, if we have a Type Error then
        # try with the compat pickle to handle subclass changes
        # pass encoding only if its not None as py2 doesn't handle
        # the param

        # cpickle
        # GH 6899
        try:
            return pkl.loads(data)
        except:
            # reg/patched pickle
            try:
                with BytesIO(data) as fh:
                    return pc.load(fh, encoding=encoding, compat=False)

            # compat pickle
            except:
                with BytesIO(data) as fh:
                    return pc.load(fh, encoding=encoding, compat=True)

    try:
        return try_read(dump)
    except:
        if PY3:
            return try_read(dump, encoding='latin1')
        raise
