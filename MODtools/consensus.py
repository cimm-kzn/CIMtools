# -*- coding: utf-8 -*-
#
# Copyright 2015 Ramil Nugmanov <stsouko@live.ru>
# This file is part of MODtools.
#
# MODtools is free software; you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Affero General Public License for more details.
#
#  You should have received a copy of the GNU Affero General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#


class ConsensusDragos(object):
    errors = dict(lad='Too few (less than %d %%) local models have applicability domains covering this structure. ',
                  diff='The other local models disagree (prediction value difference = %.2f) with the prediction of '
                       'the minority containing structure inside their applicability domain. ',
                  stp='Individual models failed to reach unanimity - prediction variance exceeds %d %% '
                      'of the property range width. ',
                  zad='None of the local models have applicability domains covering this structure. ')

    trust_desc = {5: 'Optimal', 4: 'Good', 3: 'Medium', 2: 'Low', 1: 'Bad'}
