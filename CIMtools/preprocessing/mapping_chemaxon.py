# -*- coding: utf-8 -*-
#
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
from CGRtools.containers import ReactionContainer
from .standardize import StandardizeChemAxon


class MappingChemAxon(StandardizeChemAxon):
    def __init__(self, workpath='.'):
        rules = '<?xml version="1.0" encoding="UTF-8"?><StandardizerConfiguration Version="0.1"><Actions>' \
                '<UnmapReaction ID="Unmap"/><MapReaction ID="Map Reaction" KeepMapping="false" ' \
                'MappingStyle="COMPLETE" MarkBonds="false"/></Actions></StandardizerConfiguration>'

        super().__init__(rules, workpath)

    _dtype = ReactionContainer


__all__ = ['MappingChemAxon']
