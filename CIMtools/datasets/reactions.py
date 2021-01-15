# -*- coding: utf-8 -*-
#
#  Copyright 2021 Assima Rakhimbekova <asima.rakhimbekova@outlook.com>
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
from numpy import array
from sklearn.datasets.base import _convert_data_dataframe
from sklearn.utils import Bunch
from CGRtools import RDFRead


def load_sn2(*, return_X_y=False, as_frame=False):
    """Load and return bimolecular nucleophilic substitution (SN2) reactions dataset (regression).

    ==============   ==============
    Samples total     4830
    Data              reactions, type: ReactionContainer
    Targets           real logK (-7.68) - 1.65, type: float
    ==============   ==============

    Read more in the article:
    Gimadiev, T.; Madzhidov, T.; Tetko, I.; Nugmanov, R.; Casciuc, I.;
    Klimchuk, O.; Bodrov, A.; Polishchuk, P.; Antipin, I.; Varnek, A.
    Bimolecular Nucleophilic Substitution Reactions: Predictive Models
    for Rate Constants and Molecular Reaction Pairs Analysis.
    J. Mol. Inf. 2018, 38, 1800104, doi:10.1002/minf.201800104.

    Parameters
    ----------
    return_X_y : bool, default=False
        If True, returns ``(data, target)`` instead of a Bunch object.

    as_frame : bool, default=False
        If True, the data is a pandas DataFrame including columns with
        appropriate dtypes (numeric). The target is
        a pandas DataFrame or Series depending on the number of target columns.
        If `return_X_y` is True, then (`data`, `target`) will be pandas
        DataFrames or Series as described below.

    Returns
    -------
    data : :class:`~sklearn.utils.Bunch`
        Dictionary-like object, with the following attributes.

        data : ndarray of shape (4830, )
            The data array.
        target : ndarray of shape (4830, )
            The regression target (logarithm of the reaction rate constant).

    (data, target) : tuple if ``return_X_y`` is True
    """
    data = []
    with RDFRead('SN2.rdf') as f:
        for r in f:
            r.canonicalize() # standardize dataset: (1) standardize functional groups;
            # (2) leads to kekula formula for benzene rings; (3) removes hydrogens; (4) aromatizes benzene rings.
            data.append(r)
    data = array(data)
    target = array([float(x.meta['logK']) for x in data])
    if as_frame:
        frame = None
        frame, data, target = _convert_data_dataframe(caller_name="load_sn2", data=data, target=target,
                                                      feature_names=['reaction', ], target_names=['target', ])
    if return_X_y:
        return data, target
    return Bunch(data=data, target=target)

def load_e2(*, return_X_y=False, as_frame=False):
    """Load and return bimolecular elimination (E2) reactions dataset (regression).

    ==============   ==============
    Samples total     1820
    Data              reactions, type: ReactionContainer
    Targets           real logK (-7.23) - 2.67, type: float
    ==============   ==============

    Read more in the article:
    Madzhidov, T.I.; Bodrov, A.V.; Gimadiev, T.R.; Nugmanov, R.I.; Antipin, I.S.;
    Varnek, A. Structure–reactivity relationship in bimolecular elimination
    reactions based on the condensed graph of a reaction.
    J. Struct. Chem. 2015, 56, 1227–1234, doi:10.1134/S002247661507001X.

    Parameters
    ----------
    return_X_y : bool, default=False
        If True, returns ``(data, target)`` instead of a Bunch object.

    as_frame : bool, default=False
        If True, the data is a pandas DataFrame including columns with
        appropriate dtypes (numeric). The target is
        a pandas DataFrame or Series depending on the number of target columns.
        If `return_X_y` is True, then (`data`, `target`) will be pandas
        DataFrames or Series as described below.

    Returns
    -------
    data : :class:`~sklearn.utils.Bunch`
        Dictionary-like object, with the following attributes.

        data : ndarray of shape (1820, )
            The data array.
        target : ndarray of shape (1820, )
            The regression target (logarithm of the reaction rate constant).

    (data, target) : tuple if ``return_X_y`` is True
    """
    data = []
    with RDFRead('E2.rdf') as f:
        for r in f:
            r.canonicalize() # standardize dataset: (1) standardize functional groups;
            # (2) leads to kekula formula for benzene rings; (3) removes hydrogens; (4) aromatizes benzene rings.
            data.append(r)
    data = array(data)
    target = array([float(x.meta['logK']) for x in data])
    if as_frame:
        frame = None
        frame, data, target = _convert_data_dataframe(caller_name="load_e2", data=data, target=target,
                                                      feature_names=['reaction', ], target_names=['target', ])
    if return_X_y:
        return data, target
    return Bunch(data=data, target=target)

def load_da(*, return_X_y=False, as_frame=False):
    """Load and return Diels-Alder reactions (DA) reactions dataset (regression).

    ==============   ==============
    Samples total     1866
    Data              reactions, type: ReactionContainer
    Targets           real logK (-8.511) - 8.568, type: float
    ==============   ==============

    Read more in the article:
    Madzhidov, T.I.; Gimadiev, T.R.; Malakhova, D.A.; Nugmanov, R.I.;
    Baskin, I.I.; Antipin, I.S.; Varnek, A. Structure-Reactivity modelling
    for Diels-alder reactions based on the condensed REACTION graph approach.
    J. Struct. Chem. 2017, 58, 685–691, doi:10.1134/S0022476617040023.

    Parameters
    ----------
    return_X_y : bool, default=False
        If True, returns ``(data, target)`` instead of a Bunch object.

    as_frame : bool, default=False
        If True, the data is a pandas DataFrame including columns with
        appropriate dtypes (numeric). The target is
a pandas DataFrame or Series depending on the number of target columns.
        If `return_X_y` is True, then (`data`, `target`) will be pandas
        DataFrames or Series as described below.

    Returns
    -------
    data : :class:`~sklearn.utils.Bunch`
        Dictionary-like object, with the following attributes.

        data : ndarray of shape (1866, )
            The data array.
        target : ndarray of shape (1866, )
            The regression target (logarithm of the reaction rate constant).

    (data, target) : tuple if ``return_X_y`` is True
    """
    data = []
    with RDFRead('DA.rdf') as f:
        for r in f:
            r.canonicalize() # standardize dataset: (1) standardize functional groups;
            # (2) leads to kekula formula for benzene rings; (3) removes hydrogens; (4) aromatizes benzene rings.
            data.append(r)
    data = array(data)
    target = array([float(x.meta['logK']) for x in data])
    if as_frame:
        frame = None
        frame, data, target = _convert_data_dataframe(caller_name="load_da", data=data, target=target,
                                                      feature_names=['reaction', ], target_names=['target', ])
    if return_X_y:
        return data, target
    return Bunch(data=data, target=target)


__all__ = ['load_sn2', 'load_e2', 'load_da']
