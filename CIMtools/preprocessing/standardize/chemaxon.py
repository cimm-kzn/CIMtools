# -*- coding: utf-8 -*-
#
#  Copyright 2018-2021 Ramil Nugmanov <nougmanoff@protonmail.com>
#  Copyright 2020 Zarina Ibragimova <zarinaIbr12@yandex.ru>
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
from CGRtools import MRVRead, MRVWrite, ReactionContainer
from io import StringIO, BytesIO
from logging import warning
from pandas import DataFrame
from pathlib import Path
from shutil import which
from ...base import CIMtoolsTransformerMixin
from ...exceptions import ConfigurationError


class StandardizeChemAxon(CIMtoolsTransformerMixin):
    def __init__(self, rules, *, _skip_errors=False):
        self.rules = rules
        self.__skip = _skip_errors
        self.__standardizer_obj = self.__standardizer(rules)

    def __new__(cls, *args, **kwargs):
        if cls.__standardizer is None:  # load only once
            from jnius_config import add_classpath
            add_classpath(*jars)
            from jnius import autoclass

            cls.__standardizer = autoclass('chemaxon.standardizer.Standardizer')
            cls.__importer = autoclass('chemaxon.formats.MolImporter')
            cls.__exporter = autoclass('chemaxon.formats.MolExporter')
        return super().__new__(cls)

    def __getstate__(self):
        return {'rules': self.rules}

    def __setstate__(self, state):
        super().__setstate__(state)
        self.__standardizer_obj = self.__standardizer(self.rules)

    def set_params(self, **params):
        if params:
            super().set_params(**params)
            self.__standardizer_obj = self.__standardizer(self.rules)
        return self

    def transform(self, x):
        x = super().transform(x)
        out = []
        for n, s in enumerate(x):
            with StringIO() as f:
                with MRVWrite(f) as w:
                    w.write(s)
                js = self.__importer.importMol(f.getvalue())

            try:
                self.__standardizer_obj.standardize(js)
            except Exception as e:
                if 'Invalid standardizer action' in e.args[0]:
                    raise ConfigurationError from e
                if self.__skip:
                    warning(f'structure ({n}): {s} not processed')
                    out.append([s])
                    continue
                raise ValueError from e
            mrv = self.__exporter.exportToFormat(js, 'mrv')
            with BytesIO(mrv.encode()) as f, MRVRead(f, remap=False) as r:
                p = r.read()
                if not p:
                    if self.__skip:
                        warning(f'structure ({n}): {s} export failed')
                        out.append([s])
                        continue
                    raise ValueError(f'structure ({n}): {s} export failed')
            out.append(p)
        return DataFrame(out, columns=['standardized'])

    __standardizer = None
    __importer = None
    __exporter = None


class MappingChemAxon(StandardizeChemAxon):
    def __init__(self, *, _skip_errors=False):
        rules = '<?xml version="1.0" encoding="UTF-8"?><StandardizerConfiguration Version="0.1"><Actions>' \
                '<UnmapReaction ID="Unmap"/><MapReaction ID="Map Reaction" KeepMapping="false" ' \
                'MappingStyle="COMPLETE" MarkBonds="false"/></Actions></StandardizerConfiguration>'
        super().__init__(rules, _skip_errors=_skip_errors)

    _dtype = ReactionContainer


std = which('standardize')
if std:
    jars = [str(x) for x in (Path(std).resolve().parents[1] / 'lib').iterdir() if x.name.endswith('.jar')]
    __all__ = ['StandardizeChemAxon', 'MappingChemAxon']
else:
    del StandardizeChemAxon, MappingChemAxon
    __all__ = []
