# -*- coding: utf-8 -*-
#
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
import tensorflow as tf
from tensorflow.keras import backend as K
from tensorflow.keras.layers import (Input, Lambda, Layer, Dense, Concatenate, BatchNormalization, Conv1D,
                                     TimeDistributed, Reshape, Embedding)
from tensorflow.keras.models import Model


def FTSwish(threshold=-0.2):
    def _FTSwish(x):
        return K.relu(x) * K.sigmoid(x) + threshold

    return Lambda(_FTSwish)


def p_relu():
    def _p_relu(x):
        return K.relu(x) + 0.001

    return Lambda(_p_relu)


def get_matrix_all():
    def _perm_dims(x):
        adj_m, atoms = x
        adj_m = K.expand_dims(adj_m, axis=-2)
        atoms = K.expand_dims(atoms, axis=-1)
        ans = adj_m * atoms
        return K.permute_dimensions(ans, pattern=(0, 3, 1, 2))

    return Lambda(_perm_dims)


def mask_by_adj():
    def _mask_by_adj(x):
        adj_m, x = x
        return x * K.expand_dims(adj_m, axis=-1)

    return Lambda(_mask_by_adj)


def get_top_k(k):
    def _get_top_k(x):
        f = K.permute_dimensions(x, pattern=(0, 1, 3, 2))
        f = tf.nn.top_k(f, k=k).values
        f = K.permute_dimensions(f, pattern=(0, 1, 3, 2))
        return f

    return Lambda(_get_top_k)


def adj_m_pad_k(k):
    def _adj_m_pad_k(adj_m):
        size = K.shape(adj_m)[-2]
        adj_m = tf.cond(size < k, lambda: tf.pad(adj_m, [[0, 0], [0, k - size], [0, k - size]]), lambda: adj_m)
        return adj_m

    return Lambda(_adj_m_pad_k)


def atoms_pad_k(k):
    def _atoms_pad_k(atoms):
        size = K.shape(atoms)[-2]
        atoms = tf.cond(size < k, lambda: tf.pad(atoms, [[0, 0], [0, k - size], [0, 0]]), lambda: atoms)
        return atoms

    return Lambda(_atoms_pad_k)


def connections_m_pad_k(k):
    def _connections_m_pad_k(connections_m):
        size = K.shape(connections_m)[-2]
        connections_m = tf.cond(size < k, lambda: tf.pad(connections_m, [[0, 0], [0, k - size], [0, k - size], [0, 0]]),
                                lambda: connections_m)
        return connections_m

    return Lambda(_connections_m_pad_k)


def mask_pad_by_adj():
    def _mask_pad_by_adj(x):
        adj_m_pad, x = x
        n = K.sum(K.sum(adj_m_pad, axis=-1), axis=0)
        n = tf.math.divide_no_nan(n, n)
        n = K.expand_dims(n, axis=0)
        n = K.expand_dims(n, axis=-1)
        return x * n

    return Lambda(_mask_pad_by_adj)


class Atom_Emb(Layer):
    def __init__(self, nodes_num, emb_size):
        super(Atom_Emb, self).__init__()
        self.emb = Embedding(nodes_num, emb_size, mask_zero=True)

    def call(self, inputs):
        # charge, atomic_number, is_radical = tf.split(inputs, 3, axis=-1)
        atomic_number, charge, is_radical = tf.split(inputs, 3, axis=-1)
        atomic_emb = K.squeeze(self.emb(atomic_number), axis=-2)
        return Concatenate(axis=-1)([charge, atomic_emb, is_radical])


def RMS():
    def _RMS(x):
        n = K.sum(x, axis=-1)
        n = tf.math.divide_no_nan(n, n)
        size = K.sum(n, axis=-1)
        return K.sqrt(K.sum(K.square(x), axis=1) / K.expand_dims(size, axis=-1))

    return Lambda(_RMS)


def GraphConv(n_atoms, n_connections, k=4, selector_size=20):
    adj_m = Input(shape=(None, None))
    atoms = Input(shape=(None, n_atoms))
    connections_m = Input(shape=(None, None, n_connections))

    adj_m_pad = adj_m_pad_k(k=k)(adj_m)
    atoms_pad = atoms_pad_k(k=k)(atoms)
    connections_m_pad = connections_m_pad_k(k=k)(connections_m)

    x = get_matrix_all()([adj_m_pad, atoms_pad])
    x = Concatenate(axis=-1)([x, connections_m_pad])

    selector = Dense(selector_size, kernel_initializer='he_normal')

    x = p_relu()(BatchNormalization()(selector(x)))
    x = mask_by_adj()([adj_m_pad, x])

    x = get_top_k(k=k)(x)

    ext_atoms = Lambda(lambda x: tf.pad(x, [[0, 0], [0, 0], [0, K.int_shape(connections_m_pad)[-1]]]))(atoms_pad)
    ext_atoms = Lambda(lambda x: K.expand_dims(x, axis=-2))(ext_atoms)
    ext_atoms = p_relu()(BatchNormalization()(selector(ext_atoms)))

    x = Concatenate(axis=-2)([ext_atoms, x])  # Nx(k+1)x(selector_size)

    if k % 2 != 0:
        x = TimeDistributed(Conv1D(200, kernel_size=2, kernel_initializer='he_normal'))(x)
        x = BatchNormalization()(x)
        x = FTSwish()(x)

    if k >= 4:
        x = TimeDistributed(Conv1D(200, kernel_size=k - 1, kernel_initializer='he_normal'))(x)
        x = BatchNormalization()(x)
        x = FTSwish()(x)

    x = TimeDistributed(Conv1D(100, kernel_size=3, kernel_initializer='he_normal'))(x)
    x = BatchNormalization()(x)
    x = FTSwish()(x)

    x = Reshape([-1, 100])(x)  # Nx100

    x = mask_pad_by_adj()([adj_m_pad, x])

    return Model(inputs=[adj_m, atoms, connections_m], outputs=x)


def GNN(nodes_num, connections_num, selector_size, top_k, depth):
    nodes = Input(shape=(None, 3))
    connections_matrix = Input(shape=(None, None))

    adj_m = Lambda(lambda x: tf.math.divide_no_nan(x, x))(connections_matrix)

    atoms_emb = Atom_Emb(nodes_num, 20)(nodes)

    connections_m_emb = Embedding(connections_num, 10, mask_zero=True)(connections_matrix)

    vectors = []
    for i in range(depth):
        if i == 0:
            vectors = GraphConv(n_atoms=22, n_connections=10, k=top_k, selector_size=selector_size)(
                    [adj_m, atoms_emb, connections_m_emb])
        else:
            tmp = GraphConv(n_atoms=i * 100 + 22, n_connections=10, k=top_k, selector_size=selector_size)(
                    [adj_m, Concatenate(axis=-1)([atoms_emb, vectors]), connections_m_emb])
            vectors = Concatenate(axis=-1)([vectors, tmp])

    concat_vectors = Dense(300, kernel_initializer='he_normal')(vectors)
    concat_vectors = BatchNormalization()(concat_vectors)
    concat_vectors = FTSwish()(concat_vectors)

    mols = RMS()(concat_vectors)

    return Model(inputs=[nodes, connections_matrix], outputs=mols)


__all__ = ['GNN']
