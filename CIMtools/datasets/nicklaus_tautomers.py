# -*- coding: utf-8 -*-
#
#  Copyright 2021 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from CGRtools import smiles
from numpy import array
from math import isnan
from pandas import read_excel
from pkg_resources import resource_stream
from sklearn.datasets._base import _convert_data_dataframe
from sklearn.utils import Bunch


def load_nicklaus_tautomers(*, return_X_y=False, as_frame=False, as_regression=False):
    """Load and return Nicklaus's tautomers dataset (Regression and Classification).

    ==================   ==============
    Samples total         5960
    Samples Regression    2824
    Data                  molecules, type: MoleculeContainer
    Targets               real ratio 0.0 - 1.0, type: float (Regression)
    Classes               5
    ==================   ==============

    Molecules has .meta attribute which returns dict with additional data:
    structure_id: row in original file
    tautomer_id: id of structure of tautomer in row
    additive.{n}: solvent name. {n} started from 1 id of solvent. in mixtures will be presented more additive keys.
        e.g. additive.2, additive.3 ...
    amount.{n}: amount of additive.
    prevalence (optional): Qualitative category of tautomer reported in the publication.
    temperature (optional) in Kelvin
    pH (optional)

    For Regression: The numeric proportion of tautomer based on its quantitative ratio and qualitative prevalence.

    For Classification: Quantitative ratio of tautomer compared to other tautomers.

    Numeric classification of qualitative prevalence:
    0: Not observed
    1: Less favored, less stable, minor, observed
    2: Equally, favored, major, in equilibrium,  preferred, similar spectra
    3: More favored, more stable, predominant, strongly favored
    4: Exclusively observed, only observed, only tautomer, identical tautomer

    Numeric classification of quantitative amount of tautomers:
    0: ratio = 0.0 - 0.0099
    1: ratio =  0.01 - 0.30
    2: ratio = 0.31 - 0.69
    3: ratio =  0.70 - 0.99
    4: ratio = 1

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

    as_regression : bool, default=False
        If True, returns regression subset instead of classes

    Returns
    -------
    data : :class:`~sklearn.utils.Bunch`
        Dictionary-like object, with the following attributes.

        data : ndarray of shape (n, )
            The data array.
        target : ndarray of shape (n, )
            The regression or classification target.
        feature_names: list
            The name of the dataset ('Tautomers').
        target_names: list
            The name of target (['ratio or category']).
        frame: DataFrame of shape (n, 2)
            Only present when `as_frame=True`. DataFrame with `data` and
            `target`.

    (data, target) : tuple if ``return_X_y`` is True
    """
    with resource_stream('CIMtools.datasets', 'data/tautomer_database_release_3a.xlsx') as s:
        data = read_excel(s, na_values='nul', true_values=['yes'], false_values=['no'])

    parsed = []
    for record, mol in data.iterrows():
        meta = {}
        n = m = 0
        c = mol['Size']

        sol = mol['Solvent']
        if mol['Solvent_Mixture']:
            for n, sol in enumerate(sol.split(','), 1):
                meta[f'additive.{n}'] = sol.strip()

            sol_prop = mol['Solvent_Proportion']
            if isinstance(sol_prop, str):
                sol_prop = [float(x) for x in sol_prop.split(':')]
                sum_prop = sum(sol_prop)
                for m, prop in enumerate(sol_prop, 1):
                    meta[f'amount.{m}'] = prop / sum_prop
                if n != m:
                    raise ValueError
            elif not isnan(sol_prop):
                raise ValueError
        else:
            meta['additive.1'] = sol
            meta['amount.1'] = 1.

        temp = mol['Temperature']
        if not isinstance(temp, str) and not isnan(temp):
            meta['temperature'] = temp
        ph = mol['pH']
        if not isinstance(ph, str) and not isnan(ph):
            meta['pH'] = ph

        tmp = []
        for n in range(1, c + 1):
            s = smiles(mol[f'SMILES_{n}'])
            s.kekule()
            s.standardize(fix_stereo=False)
            if s.check_valence():  # skip invalid records
                break
            s.implicify_hydrogens(fix_stereo=False)
            s.thiele()

            s.meta['structure_id'] = record + 1
            s.meta['tautomer_id'] = n
            s.meta['category'] = mol[f'Prevalence_Category_{n}']

            rat = mol[f'Quantitative_ratio_{n}']
            if not isinstance(rat, str) and not isnan(rat):
                s.meta['ratio'] = rat
            prev = mol[f'Qualitative_prevalence_{n}']
            if isinstance(prev, str):
                s.meta['prevalence'] = prev.lower()
            s.meta.update(meta)
            tmp.append(s)
        else:
            parsed.extend(tmp)

    if as_regression:
        data = array([x for x in parsed if 'ratio' in x.meta])
        target = array([x.meta['ratio'] for x in data])
        target_names = ['ratio']
    else:
        data = array(parsed)
        target = array([x.meta['category'] for x in data], dtype=int)
        target_names = ['category']

    feature_names = ['Tautomers']

    if as_frame:
        frame, data, target = _convert_data_dataframe(caller_name='load_nicklaus_tautomers', data=data, target=target,
                                                      feature_names=feature_names, target_names=target_names)
    else:
        frame = None

    if return_X_y:
        return data, target
    return Bunch(data=data, target=target, frame=frame, target_names=target_names, feature_names=feature_names)


__all__ = ['load_nicklaus_tautomers']
