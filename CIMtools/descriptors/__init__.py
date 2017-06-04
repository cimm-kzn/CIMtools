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
from collections import defaultdict
from functools import reduce
from operator import and_
from pandas import Series, concat, merge
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

    def get(self, structures, **kwargs):
        """
        :param structures: list of CGRtools data
        :param kwargs: generators specific arguments
        :return: dict(X=DataFrame, AD=Series, Y=Series, BOX=Series, structures=DataFrame)
        """
        res = defaultdict(list)

        def merge_wrap(x, y):
            return merge(x, y, how='outer', left_index=True, right_index=True)

        for gen, ad in self.__generators:
            for k, v in gen.get(structures, **kwargs).items():
                res[k].append(v)
            res['BOX'].append(Series(ad, index=res['X'][-1].columns))

        res['X'] = reduce(merge_wrap, res['X'])
        res['AD'] = reduce(and_, sorted(res['AD'], key=lambda x: len(x.index), reverse=True))
        # на данный момент не придумано как поступать с мультицентровостью. пока свойство просто дублируется.
        res['Y'] = sorted(res['Y'], key=lambda x: len(x.index), reverse=True)[0]

        res['BOX'] = concat(res['BOX'])

        if 'structures' in res:
            res['structures'] = reduce(merge_wrap, res['structures'])

        return dict(res)
