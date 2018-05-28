# -*- coding: utf-8 -*-
#
#  Copyright 2016-2018 Ramil Nugmanov <stsouko@live.ru>
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
from CGRtools.files import SDFwrite
from io import StringIO
from os import environ
from pandas import DataFrame, Series
from sklearn.base import BaseEstimator
from subprocess import run, PIPE
from .common import TransformerMixin
from ..config import UTILS_DIR, JCHEM_DIR
from ..exceptions import ConfigurationError


class Eed(BaseEstimator, TransformerMixin):
    def transform(self, x, return_domain=False):
        x = super().transform(x)

        with StringIO() as f, SDFwrite(f, mark_to_map=True) as w:
            for s in x:
                w.write(s)

            tmp = f.getvalue().encode()
        try:
            env = environ.copy()
            env.update(CLASSPATH='{}/lib/jchem.jar:{}'.format(JCHEM_DIR, UTILS_DIR))

            p = run(['java', 'Utils.react_desc', '-svm'], input=tmp, stdout=PIPE, stderr=PIPE, env=env)
        except FileNotFoundError as e:
            raise ConfigurationError(e)

        if p.returncode != 0:
            raise ConfigurationError(p.stderr.decode())

        # todo: implement error structure handler
        x, d = self.__parse_eed_output(p.stdout.decode())

        if return_domain:
            return x, d
        return x

    @classmethod
    def __parse_eed_output(cls, output):
        vector, ad = [], []
        for frag in StringIO(output):
            _, *x = frag.split()
            tmp = {}  # X vector
            for i in x:
                k, v = i.split(':')
                tmp['eed.%s' % k.rstrip()] = float(v.replace(',', '.'))
            ad.append(len(tmp) > 2)
            vector.append(tmp)

        return DataFrame(vector, columns=cls.__columns).fillna(0), Series(ad)

    __columns = ['eed.%d' % x for x in range(1, 705)]
    _dtype = MoleculeContainer
