# -*- coding: utf-8 -*-
#
#  Copyright 2020 Ramil Nugmanov <nougmanoff@protonmail.com>
#  Copyright 2020 Daniyar Mazitov <daniyarttt@gmail.com>
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
from os.path import dirname, join
from sys import modules
from ..graph_to_matrix import MoleculesToMatrix
from ...base import CIMtoolsTransformerMixin


class GNNFingerprint(CIMtoolsTransformerMixin):
    def __init__(self):
        """
        Molecules encoder
        """
        self.__m2m = MoleculesToMatrix(is_radical=True)

    def __getstate__(self):
        return {'_GNNFingerprint__m2m': self.__m2m}

    def transform(self, x):
        x = self.__m2m.transform(x).data
        x = self.__encoder(x).numpy()
        return x

    def __new__(cls, *args, **kwargs):
        if cls.__encoder is None:  # load only once
            import tensorflow as tf
            import tensorflow.keras.backend as K
            from tensorflow.keras.layers import Input, Dense
            from tensorflow.keras.models import Model
            from .gnn import GNN

            gpus = tf.config.experimental.list_physical_devices('GPU')
            if gpus:
                for gpu in gpus:
                    tf.config.experimental.set_memory_growth(gpu, True)

            atoms = Input(shape=(None, 3))
            connections_m = Input(shape=(None, None))

            m = GNN(nodes_num=119, connections_num=5, selector_size=25, top_k=4, depth=2)([atoms, connections_m])
            m = Dense(50, activation=lambda x: K.l2_normalize(x, axis=-1), kernel_initializer='truncated_normal')(m)

            encoder = Model(inputs=[atoms, connections_m], outputs=m)
            path = join(dirname(modules[__package__].__file__), 'weights.h5')
            encoder.load_weights(path)
            cls.__encoder = encoder

        return super().__new__(cls)

    __encoder = None


__all__ = ['GNNFingerprint']
