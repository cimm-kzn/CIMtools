# -*- coding: utf-8 -*-
#
#  Copyright 2019 Ramil Nugmanov <stsouko@live.ru>
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
from CGRtools.containers import MoleculeContainer, CGRContainer
from numpy import zeros, bool8
from hashlib import md5
from sklearn.base import BaseEstimator, TransformerMixin
from ..utils import iter2array


class FragmentorFingerprint(BaseEstimator, TransformerMixin):
    def __init__(self, fingerprint_size=12, bits_count=4, bits_active=2,  fragment_type=3, min_length=2, max_length=10,
                 cgr_dynbonds=0, doallways=False, useformalcharge=False, workpath='.', version='2017', verbose=False):
        """
        ISIDA Fragmentor fragments to fingerprints

        :param fingerprint_size: exponent of 2 of fingerprint length
        :param bits_count: include number of fragment descriptors into fingerprint. for example by default:
            if number is one: only one set of bits will be activated.
            if number is tho: additional set of bits will be activated (totally up to 2*bits_active).
            if number is four or greater: will be activated up 4*bits_active bits
        :param bits_active: number of activated bits for each fragment (need for prevent collision bit lost)
        :param workpath: path for temp files.
        :param version: fragmentor version. need for selecting Fragmentor executables named as fragmentor-{version}
        """
        self.__fragmentor = Fragmentor(fragment_type=fragment_type, min_length=min_length, max_length=max_length,
                                       cgr_dynbonds=cgr_dynbonds, doallways=doallways, useformalcharge=useformalcharge,
                                       header=False, workpath=workpath, version=version, verbose=verbose)

        self.fingerprint_size = fingerprint_size
        self.bits_count = bits_count
        self.bits_active = bits_active
        self.fragment_type = fragment_type
        self.min_length = min_length
        self.max_length = max_length
        self.cgr_dynbonds = cgr_dynbonds
        self.doallways = doallways
        self.useformalcharge = useformalcharge
        self.workpath = workpath
        self.version = version
        self.verbose = verbose

    def __getstate__(self):
        return {k: v for k, v in super().__getstate__().items()
                if k not in ('_FragmentorFingerprint__fragmentor', 'workpath')}

    def __setstate__(self, state):
        super().__setstate__(state)
        self.__fragmentor = Fragmentor(**{k: v for k, v in state.items()
                                          if k not in ('fingerprint_size', 'bits_count', 'bits_active')})
        self.workpath = '.'

    def __del__(self):
        self.__fragmentor.delete_work_path()

    def set_params(self, **params):
        if not params:
            return self

        super().set_params(**params)
        self.__fragmentor.set_params(**{k: v for k, v in params.items()
                                        if k not in ('fingerprint_size', 'bits_count', 'bits_active')})
        return self

    def set_work_path(self, workpath):
        self.__fragmentor.set_work_path(workpath)

    def delete_work_path(self):
        self.__fragmentor.delete_work_path()

    def fit(self, x, y=None):
        return self

    def transform(self, x):
        x = self.transform_bitset(x)
        out = zeros((len(x), 2 ** self.fingerprint_size), dtype=bool8)
        for i, j in enumerate(x):
            out[i, j] = True
        return out

    def transform_bitset(self, x):
        x = iter2array(x, dtype=(MoleculeContainer, CGRContainer))
        mask = 2 ** self.fingerprint_size - 1
        fp_count = self.bits_count
        fp_active = self.bits_active * 2

        df = self.__fragmentor.transform(x)
        bits_map = {}
        for f in df.columns:
            prev = []
            for i in range(1, fp_count + 1):
                bs = md5(f'{i}_{f}'.encode()).digest()
                bits_map[(f, i)] = prev = [int.from_bytes(bs[r: r + 2], 'big') & mask
                                           for r in range(0, fp_active, 2)] + prev

        out = []
        for s, (_, row) in zip(x, df.iterrows()):
            active_bits = set()
            # add atomic bits
            if isinstance(s, MoleculeContainer):
                for i in set(int(a) for _, a in s.atoms()):
                    active_bits.add(i & mask)
                    active_bits.add((i >> 5) & mask)  # charge and isotope excluded
            else:
                for i in set(int(a) for _, a in s.atoms()):
                    active_bits.add(i & mask)
                    active_bits.add((i >> 10) & mask)  # charge and isotope excluded
            # add fragment bits
            for k, v in zip(df.columns, row):
                if v > fp_count:
                    active_bits.update(bits_map[(k, fp_count)])
                elif v:
                    active_bits.update(bits_map[(k, v)])
            out.append(list(active_bits))
        return out


try:
    from .fragmentor import Fragmentor
    __all__ = ['FragmentorFingerprint']
except ImportError:
    del FragmentorFingerprint
    __all__ = []
