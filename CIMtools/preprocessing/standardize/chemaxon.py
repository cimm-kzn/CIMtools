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
from CGRtools.files import MRVread, MRVwrite
from io import StringIO, BytesIO
from requests import post
from requests.exceptions import RequestException
from subprocess import run, PIPE
from sklearn.base import BaseEstimator
from ..common import iter2array, TransformerMixin
from ...config import STANDARDIZER, CHEMAXON
from ...exceptions import ConfigurationError


class StandardizeChemAxon(BaseEstimator, TransformerMixin):
    def __init__(self, rules):
        self.rules = rules

    def transform(self, x, y=None):
        x = super().transform(x, y)
        return iter2array(self.__processor_m(x) if x.size > 1 else self.__processor_s(x[0]), allow_none=True)

    def __processor_m(self, structures):
        with StringIO() as f, MRVwrite(f) as w:
            for s in structures:
                w.write(s)
            w.finalize()
            tmp = f.getvalue().encode()
        try:
            p = run([STANDARDIZER, '-c', self.rules, '-f', 'mrv', '-g'], input=tmp, stdout=PIPE, stderr=PIPE)
        except FileNotFoundError as e:
            raise ConfigurationError(e)

        if p.returncode != 0:
            raise ConfigurationError(p.stderr.decode())

        with BytesIO(p.stdout) as f, MRVread(f) as r:
            res = r.read()
            for x in p.stderr.decode().split('\n'):
                if x.startswith('SEVERE: Error at molecule No.'):
                    res.insert(int(x[29:].split(':', 1)[0]) - 1, None)
            return res

    def __processor_s(self, structure):
        with StringIO() as f, MRVwrite(f) as w:
            w.write(structure)
            w.finalize()
            data = dict(structure=f.getvalue(), parameters='mrv',
                        filterChain=[dict(filter='standardizer',
                                          parameters=dict(standardizerDefinition=self.rules))])
        try:
            res = self.__chemaxon_rest(data)
        except RequestException as e:
            raise ConfigurationError(e)

        if not res:
            return [None]

        with BytesIO(res['structure'].encode()) as f, MRVread(f) as r:
            return r.read()

    @classmethod
    def __chemaxon_rest(cls, data):
        q = post(cls.__rest_url, json=data, headers={'content-type': 'application/json'}, timeout=20)
        if q.status_code in (201, 200):
            return q.json()

    __rest_url = '%s/rest-v0/util/calculate/molExport' % CHEMAXON
