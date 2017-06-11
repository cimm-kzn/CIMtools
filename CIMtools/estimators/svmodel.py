# -*- coding: utf-8 -*-
#
#  Copyright 2015-2017 Ramil Nugmanov <stsouko@live.ru>
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
from functools import partial
from itertools import product
from sklearn.svm import SVR, SVC
from .basemodeler import BaseModel, MinMaxScalerWrapper, FitContainer


class SVModel(BaseModel):
    def __init__(self, descriptors_generator, structures, fit_params=None, estimator='svr', probability=False,
                 nfold=5, repetitions=1, dispcoef=0, fit='rmse', scorers=('rmse', 'r2'), normalize=False,
                 domain=None, domain_params=None, domain_normalize=False, n_jobs=2, max_iter=100000, **kwargs):
        self.__max_iter = max_iter
        self.__estimator = estimator
        self.__probability = [probability]

        BaseModel.__init__(self, descriptors_generator, structures, fit_params=fit_params, repetitions=repetitions,
                           nfold=nfold, dispcoef=dispcoef, fit=fit, scorers=scorers, normalize=normalize,
                           domain=domain, domain_params=domain_params, domain_normalize=domain_normalize,
                           n_jobs=n_jobs, **kwargs)

    __estimators = dict(svr=SVR, svc=SVC)

    def _init_unpickle(self, models, scalers, scores, fitparams, pred, prob, **config):
        super(SVModel, self)._init_unpickle(**config)
        self._model = FitContainer(models=[x for x in models],
                                   scalers=[(x and MinMaxScalerWrapper.unpickle(x) or None) for x in scalers],
                                   params=fitparams, pred=pred, prob=prob, scores=scores)

    def pickle(self):
        config = super(SVModel, self).pickle()
        config.update(models=[x for x in self._model.models])
        return config

    @classmethod
    def unpickle(cls, config):
        if 'models' not in config:
            raise Exception('Invalid config')
        BaseModel.unpickle(config)
        obj = cls.__new__(cls)
        obj._init_unpickle(**config)
        return obj

    @property
    def _estimator(self):
        return partial(self.__estimators[self.__estimator], max_iter=self.__max_iter)

    def _prepare_params(self, param):
        base = dict(C=param['C'], tol=param['tol'])
        base.update(dict(epsilon=param['epsilon'])
                    if self.__estimator == 'svr' else dict(probability=self.__probability))

        if param['kernel'][0] == 'linear':  # u'*v
            base.update(kernel=['linear'])
        elif param['kernel'][0] == 'rbf':  # exp(-gamma*|u-v|^2)
            base.update(kernel=['rbf'], gamma=param['gamma'])
        elif param['kernel'][0] == 'sigmoid':  # tanh(gamma*u'*v + coef0)
            base.update(kernel=['sigmoid'], gamma=param['gamma'], coef0=param['coef0'])
        elif param['kernel'][0] == 'poly':  # (gamma*u'*v + coef0)^degree
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
