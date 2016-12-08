# -*- coding: utf-8 -*-
#
# Copyright 2015 Ramil Nugmanov <stsouko@live.ru>
# This file is part of PREDICTOR.
#
# PREDICTOR is free software; you can redistribute it and/or modify
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
import json
import requests
from .config import CHEMAXON, PREDICTOR
from MWUI.config import AdditiveType


def serverget(url, params):
    for _ in range(2):
        try:
            q = requests.get("%s/api/%s" % (PREDICTOR, url), params=params, timeout=20)
        except:
            continue
        else:
            if q.status_code in (201, 200):
                return q.json()
            else:
                continue
    else:
        return False


def serverpost(url, data):
    for _ in range(2):
        try:
            q = requests.post("%s/api/%s" % (PREDICTOR, url), json=data, timeout=20)
        except:
            continue
        else:
            if q.status_code in (201, 200):
                return q.json()
            else:
                continue
    else:
        return False


def chemaxpost(url, data):
    for _ in range(2):
        try:
            q = requests.post("%s/rest-v0/util/%s" % (CHEMAXON, url), data=json.dumps(data),
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


def get_additives():
    res = serverget('resources/additives', None)
    if res:
        for a in res:
            a['type'] = AdditiveType(a['type'])
    return res or []
