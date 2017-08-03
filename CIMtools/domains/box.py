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
from .basedomain import Domain
from ..pandas import to_pickle, read_pickle


class Box(Domain):
    def __init__(self, *_):
        """
        Box don't has fit params.
        """
        pass

    def pickle(self):
        return dict(min=to_pickle(self.__x_min), max=to_pickle(self.__x_max))

    @classmethod
    def unpickle(cls, config):
        if {'min', 'max'}.difference(config):
            raise Exception('Invalid config')
        obj = cls.__new__(cls)
        obj._init_unpickle(read_pickle(config['min']), read_pickle(config['max']))
        return obj

    def _init_unpickle(self, _min, _max):
        self.__x_min = _min
        self.__x_max = _max

    def fit(self, x, y=None):
        self.__x_min = x.min()
        self.__x_max = x.max()

    def predict(self, x):
        return ((x - self.__x_min).min(axis=1) >= 0) & ((self.__x_max - x).min(axis=1) >= 0)

    @staticmethod
    def _prepare_params(_):
        pass

    __x_min = __x_max = None
