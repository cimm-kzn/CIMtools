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
        self.__init_common(estimator, max_iter, probability)
        BaseModel.__init__(self, descriptors_generator, structures, fit_params=fit_params, repetitions=repetitions,
                           nfold=nfold, dispcoef=dispcoef, fit=fit, scorers=scorers, normalize=normalize,
                           domain=domain, domain_params=domain_params, domain_normalize=domain_normalize,
                           n_jobs=n_jobs, **kwargs)

    def _init_unpickle(self, models, estimator, max_iter, probability, scalers, scores, fitparams, pred, prob,
                       **config):
        self.__init_common(estimator, max_iter, probability)
        super(SVModel, self)._init_unpickle(**config)

        svm = []
        for x in models:
            m = self._estimator(**fitparams)
            for k, v in x.items():
                setattr(m, k, v)
            svm.append(m)

        self._model = FitContainer(models=svm, params=fitparams, pred=pred, prob=prob, scores=scores,
                                   scalers=[(x and MinMaxScalerWrapper.unpickle(x) or None) for x in scalers])

    def __init_common(self, estimator, max_iter, probability):
        self.__max_iter = max_iter
        self.__estimator = estimator
        self.__probability = probability

    def pickle(self):
        config = super(SVModel, self).pickle()
        config.update(estimator=self.__estimator, max_iter=self.__max_iter, probability=self.__probability,
                      models=[dict(shape_fit_=x.shape_fit_, support_=x.support_, support_vectors_=x.support_vectors_,
                                   n_support_=x.n_support_, _dual_coef_=x._dual_coef_, _intercept_=x._intercept_,
                                   probA_=x.probA_, probB_=x.probB_, _sparse=x._sparse, _gamma=x._gamma)
                              for x in self._model.models])
        return config

    @classmethod
    def unpickle(cls, config):
        if {'models', 'estimator', 'max_iter', 'probability'}.difference(config):
            raise Exception('Invalid config')
        BaseModel.unpickle(config)
        obj = cls.__new__(cls)
        obj._init_unpickle(**config)
        return obj

    @property
    def _estimator(self):
        if self.__estimator == 'svc':
            return partial(self.__estimators[self.__estimator], max_iter=self.__max_iter,
                           probability=self.__probability)
        return partial(self.__estimators[self.__estimator], max_iter=self.__max_iter)

    def _prepare_params(self, param):
        base = dict(C=param['C'], tol=param['tol'])
        if self.__estimator == 'svr':
            base['epsilon'] = param['epsilon']

        if param['kernel'][0] == 'linear':  # u'*v
            base['kernel'] = ['linear']
        elif param['kernel'][0] == 'rbf':  # exp(-gamma*|u-v|^2)
            base.update(kernel=['rbf'], gamma=param['gamma'])
        elif param['kernel'][0] == 'sigmoid':  # tanh(gamma*u'*v + coef0)
            base.update(kernel=['sigmoid'], gamma=param['gamma'], coef0=param['coef0'])
        elif param['kernel'][0] == 'poly':  # (gamma*u'*v + coef0)^degree
            base.update(kernel=['poly'], gamma=param['gamma'], coef0=param['coef0'], degree=param['degree'])

        elif isinstance(param['kernel'], list) and all(callable(x) for x in param['kernel']):
            base['kernel'] = param['kernel']
        elif callable(param['kernel']):
            base['kernel'] = [param['kernel']]
        else:
            return []

        k_list = []
        v_list = []
        for k, v in base.items():
            k_list.append(k)
            v_list.append(v)

        return [{k: v for k, v in zip(k_list, x)} for x in product(*v_list)]

    __estimators = dict(svr=SVR, svc=SVC)
