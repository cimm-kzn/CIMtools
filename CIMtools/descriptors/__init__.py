# -*- coding: utf-8 -*-
#
#  Copyright 2016, 2017 Ramil Nugmanov <stsouko@live.ru>
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
from functools import reduce
from operator import and_
from pandas import merge
from .basegenerator import DataContainer
from .descriptorsdict import DescriptorsDict
from .eed import Eed
from .fragmentor import Fragmentor
from .cxcalc import Pkab


class DescriptorsChain(object):
    def __init__(self, *args):
        """
        chaining multiple descriptor generators.
        concatenate X vectors and merge AD
        :param args: set of generators
        """
        self.__generators = args

    def pickle(self):
        return [dict(name=gen.__class__.__name__, config=gen.pickle()) for gen in self.__generators]

    @classmethod
    def unpickle(cls, config):
        if not isinstance(config, list):
            raise Exception('Invalid config')
        return DescriptorsChain(*(cls._generators[x['name']].unpickle(x['config']) for x in config))

    def set_work_path(self, workpath):
        for gen in self.__generators:
            if hasattr(gen, 'set_work_path'):
                gen.set_work_path(workpath)

    def delete_work_path(self):
        for gen in self.__generators:
            if hasattr(gen, 'delete_work_path'):
                gen.delete_work_path()

    def get(self, structures, **kwargs):
        """
        :param structures: list of CGRtools data
        :param kwargs: generators specific arguments
        :return: DataContainer
        """
        _x, _y, _ad = [], [], []

        for gen in self.__generators:
            out = gen.get(structures, **kwargs)
            _x.append(out.X)
            _y.append(out.Y)
            _ad.append(out.AD)

        _x = reduce(self.__merge_wrap, _x)
        _ad = reduce(and_, sorted(_ad, key=lambda x: len(x.index), reverse=True))
        # на данный момент не придумано как поступать с мультицентровостью. пока свойство просто дублируется.
        _y = sorted(_y, key=lambda x: len(x.index), reverse=True)[0]

        return DataContainer(_x, _y, _ad)

    @staticmethod
    def __merge_wrap(x, y):
        return merge(x, y, how='outer', left_index=True, right_index=True)

    _generators = {DescriptorsDict.__name__: DescriptorsDict, Fragmentor.__name__: Fragmentor,
                   Pkab.__name__: Pkab, Eed.__name__: Eed}


__all__ = [DescriptorsChain.__name__] + list(DescriptorsChain._generators)
