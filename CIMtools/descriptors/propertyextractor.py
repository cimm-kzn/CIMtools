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


class PropertyExtractor(object):
    def __init__(self, name):
        self.__name = name

    def get_property(self, meta, marks=None):
        """
        for marked atom property can named property_name.1-2-3 - where 1-2-3 sorted marked atoms.
        for correct work NEED in rdf mapping started from 1 without breaks.
        or used common property with key property_name
        :param meta: dict of data
        :param marks: list of marked atoms
        :return: float property or None
        """
        tmp = marks and meta.get('%s.%s' % (self.__name,
                                            '-'.join(str(x) for x in sorted(marks)))) or meta.get(self.__name)

        return float(tmp) if tmp else None
