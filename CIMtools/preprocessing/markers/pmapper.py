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
from CGRtools.containers import MoleculeContainer
from CGRtools.files import RDFread, RDFwrite, SDFwrite
from io import StringIO
from itertools import product
from os import close
from pathlib import Path
from subprocess import run, PIPE, STDOUT
from tempfile import mkstemp
from xml.etree import ElementTree
from sklearn.base import BaseEstimator
from ..common import iter2array, nested_iter_to_2d_array, TransformerMixin
from ..config import PMAPPER
from ...exceptions import ConfigurationError


def remove_namespace(doc, namespace):
    """Remove namespace in the passed document in place."""
    ns = u'{%s}' % namespace
    nsl = len(ns)
    for elem in doc.getiterator():
        if elem.tag.startswith(ns):
            elem.tag = elem.tag[nsl:]
    return doc


class AtomMarkerPharmacophore(BaseEstimator, TransformerMixin):
    def __init__(self, marker_rules, workpath='.'):
        self.marker_rules = marker_rules
        self.set_work_path(workpath)
        self.__load()

    def __getstate__(self):
        return {k: v for k, v in super().__getstate__().items() if not k.startswith('_AtomMarkerPharmacophore__')}

    def __setstate__(self, state):
        super().__setstate__(state)
        self.set_work_path('.')
        self.__load()

    def __load(self):
        self.__markers = list(set(x.get('Symbol') for x in remove_namespace(ElementTree.fromstring(self.marker_rules),
                                                                            'http://www.chemaxon.com').iter('AtomSet')))

    def set_work_path(self, workpath):
        self.delete_work_path()

        fd, fn = mkstemp(prefix='clr_', suffix='.xml', dir=workpath)
        self.__config = Path(fn)
        with self.__config.open('w') as f:
            f.write(self.marker_rules)
        close(fd)

    def delete_work_path(self):
        if self.__config is not None:
            self.__config.unlink()
            self.__config = None

    def get_count(self):
        return len(self.__markers)

    def transform(self, x):
        """
        marks atoms in 7th col of sdf.
        if molecule has atom mapping - will be used mapping.
        """
        x = super().transform(x)

        with StringIO() as f, SDFwrite(f) as w:
            for s in x:
                w.write(s)
            tmp = f.getvalue().encode()

        try:
            p = run([PMAPPER, '-c', str(self.__config)], input=tmp, stdout=PIPE, stderr=PIPE)
        except FileNotFoundError as e:
            raise ConfigurationError(e)

        if p.returncode != 0:
            raise ConfigurationError(p.stderr.decode())

        marks = p.stdout.decode().split()
        output = []
        for s, mark in zip(x, marks):
            found = [[] for _ in range(self.get_count())]
            for n, m in zip(s.nodes(), mark.split(';')):
                if m:
                    tmp = s.copy()
                    tmp.node[n]['mark'] = '1'
                    found[self.__markers.index(m)].append((n, tmp))

            for x in found:
                if not x:
                    x.append((None, s))

            output.append([list(x) for x in product(*found)])
        return output

    __config = None
    _dtype = MoleculeContainer
