# -*- coding: utf-8 -*-
#
# Copyright 2015, 2016 Ramil Nugmanov <stsouko@live.ru>
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
from itertools import product
from sklearn.svm import SVR, SVC
from .basemodeler import BaseModel


class SVModel(BaseModel):
    def __init__(self, descriptorgen, fitparams, structures, workpath='.', nfold=5, repetitions=1, rep_boost=100,
                 dispcoef=0, fit='rmse', estimator='svr', probability=False, scorers=('rmse', 'r2'),
                 normalize=False, n_jobs=2, **kwargs):

        self.estimator = self.__estimators[estimator]
        self.__estimator = estimator
        self.__probability = [probability]
        self.fitparams = fitparams
        BaseModel.__init__(self, descriptorgen, structures, nfold=nfold, repetitions=repetitions,
                           rep_boost=rep_boost, dispcoef=dispcoef, fit=fit, scorers=scorers, workpath=workpath,
                           normalize=normalize, n_jobs=n_jobs, **kwargs)

    __estimators = dict(svr=SVR, svc=SVC)

    def prepareparams(self, param):
        base = dict(C=param['C'], tol=param['tol'])
        base.update(dict(epsilon=param['epsilon'])
                    if self.__estimator == 'svr' else dict(probability=self.__probability))

        if param['kernel'] == 'linear':  # u'*v
            base.update(kernel=['linear'])
        elif param['kernel'] == 'rbf':  # exp(-gamma*|u-v|^2)
            base.update(kernel=['rbf'], gamma=param['gamma'])
        elif param['kernel'] == 'sigmoid':  # tanh(gamma*u'*v + coef0)
            base.update(kernel=['sigmoid'], gamma=param['gamma'], coef0=param['coef0'])
        elif param['kernel'] == 'poly':  # (gamma*u'*v + coef0)^degree
            base.update(kernel=['poly'], gamma=param['gamma'], coef0=param['coef0'], degree=param['degree'])

        elif isinstance(param['kernel'], list) and all(callable(x) for x in param['kernel']):
            base.update(kernel=param['kernel'])
        elif callable(param['kernel']):
            base.update(kernel=[param['kernel']])
        else:
            return []

        k_list = []
        v_list = []
        for k, v in base.items():
            k_list.append(k)
            v_list.append(v)

        return [{k: v for k, v in zip(k_list, x)} for x in product(*v_list)]
