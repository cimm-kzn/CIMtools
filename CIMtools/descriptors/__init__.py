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
from pandas import Series, concat, merge
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
        :param args: set of generators or set of list[generator, consider their AD {True|False}]
        """
        if isinstance(args[0], tuple):
            self.__generators = args
        else:
            self.__generators = [(x, True) for x in args]

    def pickle(self):
        tmp = {}
        for gen, box in self.__generators:
            tmp[gen.__class__.__name__] = [gen.pickle(), box]

        return tmp

    @classmethod
    def unpickle(cls, config):
        return DescriptorsChain(*((globals()[k].unpickle(gc), b) for k, (gc, b) in config.items()))

    def set_work_path(self, workpath):
        for gen, _ in self.__generators:
            if hasattr(gen, 'set_work_path'):
                gen.set_work_path(workpath)

    def delete_work_path(self):
        for gen, _ in self.__generators:
            if hasattr(gen, 'delete_work_path'):
                gen.delete_work_path()

    def get(self, structures, return_box=False, **kwargs):
        """
        :param structures: list of CGRtools data
        :param return_box: return bitmap of descriptors for box AD  
        :param kwargs: generators specific arguments
        :return: dict(X=DataFrame, AD=Series, Y=Series, BOX=Series, structures=DataFrame)
        """
        _x, _y, _ad, _box = [], [], [], []

        for gen, ad in self.__generators:
            out = gen.get(structures, **kwargs)
            _x.append(out.X)
            _y.append(out.Y)
            _ad.append(out.AD)
            _box.append(Series(ad, index=out.X.columns))

        _x = reduce(self.__merge_wrap, _x)
        _ad = reduce(and_, sorted(_ad, key=lambda x: len(x.index), reverse=True))
        # на данный момент не придумано как поступать с мультицентровостью. пока свойство просто дублируется.
        _y = sorted(_y, key=lambda x: len(x.index), reverse=True)[0]

        _box = concat(_box)
        out = DataContainer(_x, _y, _ad)
        return (out, _box) if return_box else out

    @staticmethod
    def __merge_wrap(x, y):
        return merge(x, y, how='outer', left_index=True, right_index=True)


__all__ = [DescriptorsChain.__name__, DescriptorsDict.__name__, Fragmentor.__name__, Pkab.__name__, Eed.__name__]
