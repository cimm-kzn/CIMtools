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
from requests import post
from ..config import CHEMAXON


def chemax_post(url, data):
    for _ in range(2):
        try:
            q = post("%s/rest-v0/util/%s" % (CHEMAXON, url), data=json.dumps(data),
                     headers={'content-type': 'application/json'}, timeout=20)
        except:
            continue
        else:
            if q.status_code in (201, 200):
                return q.json()
            else:
                continue
    else:
        return False
