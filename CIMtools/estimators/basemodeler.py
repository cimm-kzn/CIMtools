# -*- coding: utf-8 -*-
#
#  Copyright 2016, 2017 Ramil Nugmanov <stsouko@live.ru>
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
from abc import ABC, abstractmethod
from collections import defaultdict, namedtuple
from itertools import count
from math import sqrt
from numpy import inf, mean, var, arange
from pandas import DataFrame, Series, concat
from sklearn.model_selection import KFold
from sklearn.externals.joblib import Parallel, delayed
from sklearn.metrics import r2_score, cohen_kappa_score
from sklearn.preprocessing import MinMaxScaler
from sklearn.utils import shuffle
from ..descriptors import *
from ..domains import *
from ..scorers import *


FoldContainer = namedtuple('FoldContainer', ['estimator', 'scaler', 'pred', 'prob'])
FitContainer = namedtuple('FitContainer', ['models', 'scalers', 'params', 'pred', 'prob', 'scores'])
ResultContainer = namedtuple('ResultContainer', ['prediction', 'probability', 'domain'])


class MinMaxScalerWrapper(MinMaxScaler):
    def pickle(self):
        return dict(min=self.min_, scale=self.scale_, copy=self.copy)

    @classmethod
    def unpickle(cls, config):
        args = {'min', 'scale', 'copy'}
        if args.difference(config):
            raise Exception('Invalid config')
        obj = cls.__new__(cls)
        obj.min_ = config['min']
        obj.copy = config['copy']
        obj.scale_ = config['scale']
        return obj


def _kfold(est, est_params, x_train, y_train, x_test, normalize):
    if normalize:
        normal = MinMaxScalerWrapper()
        x_train = DataFrame(normal.fit_transform(x_train), columns=x_train.columns, index=x_train.index)
        x_test = DataFrame(normal.transform(x_test), columns=x_test.columns, index=x_test.index)
    else:
        normal = None

    model = est(**est_params)
    model.fit(x_train, y_train)

    y_pred = Series(model.predict(x_test), index=x_test.index)
    y_prob = DataFrame(model.predict_proba(x_test),
                       index=x_test.index, columns=model.classes_) if hasattr(model, 'predict_proba') else None

    return FoldContainer(estimator=model, scaler=normal, pred=y_pred, prob=y_prob)


class BaseModel(ABC):
    def __init__(self, descriptors_generator, structures, fit_params=None, nfold=5, repetitions=1, dispcoef=0,
                 fit='rmse', scorers=('rmse', 'r2'), normalize=False,
                 domain=None, domain_params=None, domain_normalize=False, n_jobs=2, **kwargs):
        self.__fit_params = fit_params
        self.__generator = descriptors_generator

        print("Descriptors generation start")
        xy = descriptors_generator.get(structures, **kwargs)
        self.__x = xy.X
        self.__y = xy.Y
        print("Descriptors generated")
        print('========================================\n'
              'Y mean +- variance = %s +- %s\n'
              '  max = %s, min = %s\n'
              '========================================' %
              (self.__y.mean(), sqrt(self.__y.var()), self.__y.max(), self.__y.min()))

        self.__domain_class = domain or Box
        self.__init_common(nfold, repetitions, normalize, domain_normalize)

        # need for fit only
        self.__domain_params = domain_params
        self.__cv = self.__cv_naive
        self.__n_jobs = n_jobs
        self.__disp_coef = dispcoef
        self.__scorers = {x: self.__scorers_dict[x] for x in scorers}
        self.__fit_score = '$' + (fit if fit in scorers else scorers[0])
        self.__score_reporter = '\n'.join(['{0} +- variance = %({0})s +- %({0}_var)s'.format(i) for i in scorers])
        self.__scores_unfitted = FitContainer(None, None, None, None, None, scores={'$%s' % x: inf for x in scorers})

        self._model, self.__domain = self.__fit()

    def _init_unpickle(self, generator, generator_class, normalize, repetitions, nfold, x, y, domain_class, domain,
                       domain_normalize, domain_scalers, domain_scores, domain_fitparams, domain_pred, domain_prob):
        self.__x = x
        self.__y = y
        self.__generator = globals()[generator_class].unpickle(generator)
        self.__domain_class = globals()[domain_class]
        self.__domain = FitContainer(models=[self.__domain_class.unpickle(x) for x in domain],
                                     scalers=[(x and MinMaxScalerWrapper.unpickle(x) or None) for x in domain_scalers],
                                     params=domain_fitparams, pred=domain_pred, prob=domain_prob, scores=domain_scores)
        self.__init_common(nfold, repetitions, normalize, domain_normalize)

    def __init_common(self, nfold, repetitions, normalize, domain_normalize):
        self.__nfold = nfold
        self.__repetitions = repetitions
        self.__normalize = normalize
        self.__domain_normalize = domain_normalize

        self.__pickle = dict(normalize=normalize, repetitions=repetitions, nfold=nfold, x=self.__x, y=self.__y,
                             domain_normalize=domain_normalize)

    @abstractmethod
    def pickle(self):
        config = self.__pickle.copy()
        config.update(generator=self.__generator.pickle(), generator_class=self.__generator.__class__.__name__,
                      domain=[x.pickle() for x in self.__domain.models], domain_class=self.__domain_class.__name__,
                      domain_scalers=[(x and x.pickle() or None) for x in self.__domain.scalers],
                      domain_scores=self.__domain.scores, domain_fitparams=self.__domain.params,
                      domain_pred=self.__domain.pred, domain_prob=self.__domain.prob,
                      scalers=[(x and x.pickle() or None) for x in self._model.scalers],
                      scores=self._model.scores, fitparams=self._model.params,
                      pred=self._model.pred, prob=self._model.prob)
        return config

    @classmethod
    @abstractmethod
    def unpickle(cls, config):
        args = {'generator', 'generator_class', 'domain', 'domain_normalize', 'domain_scalers', 'domain_scores',
                'scalers', 'scores', 'normalize', 'repetitions', 'nfold', 'x', 'y', 'models', 'domain_fitparams',
                'fitparams', 'domain_pred', 'domain_prob', 'pred', 'prob', 'domain_class'}
        if args.difference(config):
            raise Exception('Invalid config')

    def predict(self, structures, **kwargs):
        res = self.__generator.get(structures, **kwargs)
        d_x, d_ad = res.X, res.AD

        pred, prob, dom = [], [], []
        for i, (model, scaler, domain, domain_scaler) in enumerate(zip(self._model.models, self._model.scalers,
                                                                       self.__domain.models, self.__domain.scalers)):
            x_t = DataFrame(scaler.transform(d_x), columns=d_x.columns) if scaler else d_x
            pred.append(Series(model.predict(x_t), index=d_x.index))

            if hasattr(model, 'predict_proba'):
                y_pa = DataFrame(model.predict_proba(x_t), index=d_x.index, columns=model.classes_)
                prob.append(y_pa)

            x_t = DataFrame(domain_scaler.transform(d_x), columns=d_x.columns) if domain_scaler else d_x
            dom.append(Series(domain.predict(x_t), index=d_x.index))

        return ResultContainer(prediction=concat(pred, axis=1), domain=concat(dom, axis=1),
                               probability=concat(prob, axis=1, keys=range(len(prob))) if prob else None)

    @abstractmethod
    def _prepare_params(self, param):
        return []

    @property
    @abstractmethod
    def _estimator(self):
        pass

    def set_work_path(self, workpath):
        if hasattr(self.__generator, 'set_work_path'):
            self.__generator.set_work_path(workpath)

    def delete_work_path(self):
        if hasattr(self.__generator, 'delete_work_path'):
            self.__generator.delete_work_path()

    def get_model_stats(self):
        stat = dict(fitparams=self._model.params, repetitions=self.__repetitions, nfolds=self.__nfold,
                    normalize=self.__normalize, dragostolerance=sqrt(self.__y.var()),
                    domain=self.__domain_class.__name__, domain_normalize=self.__domain_normalize)

        stat.update({x: y for x, y in self._model.scores.items() if not x.startswith('$')})
        return stat

    def get_fit_predictions(self):
        return dict(property=self.__y, prediction=self._model.pred, probability=self._model.prob,
                    domain=self.__domain.pred)

    def get_descriptors_generator(self):
        return self.__generator

    def __fit(self):
        depindex, maxdep, fitparams, max_params = self.__inception(self.__fit_params)
        fcount = count(0)
        bestmodel = badmodel = self.__scores_unfitted
        for param, md, di in zip(fitparams, maxdep, depindex):
            var_kern_model = badmodel
            while param:
                var_param_model = badmodel
                tmp = self._prepare_params({k: list(v) for k, v in param.items()})
                for i in tmp:
                    print('%d (max=%d): fit model with params:' % (next(fcount), max_params), i)
                    fitted = self.__cv_model_fit(i)
                    print(self.__score_reporter % fitted.scores)
                    if fitted.scores[self.__fit_score] < var_param_model.scores[self.__fit_score]:
                        var_param_model = fitted

                if var_param_model.scores[self.__fit_score] < var_kern_model.scores[self.__fit_score]:
                    var_kern_model = var_param_model
                    param = self.__dive(var_kern_model.params, param, md, di)
                else:
                    break
            if var_kern_model.scores[self.__fit_score] < bestmodel.scores[self.__fit_score]:
                bestmodel = var_kern_model

        print('========================================\nparams:', bestmodel.params)
        print(self.__score_reporter % bestmodel.scores)
        print('%s variants checked\n========================================' % next(fcount))

        if self.__domain_params:
            raise Exception('NOT IMPLEMENTED')
        else:
            bestdomain = self.__domain_fit()

        return bestmodel, bestdomain

    def __cv_model_fit(self, fitparams):
        models, scalers, y_pred, y_prob = [], [], [], []
        fold_scorers = defaultdict(list)
        parallel = Parallel(n_jobs=self.__n_jobs)
        folds = parallel(delayed(_kfold)(self._estimator, fitparams, self.__x.iloc[train], self.__y.iloc[train],
                                         self.__x.iloc[test], self.__normalize) for train, test in self.__cv())

        #  street magic. split folds to repetitions
        for kfold in zip(*[iter(folds)] * self.__nfold):
            ky_pred, ky_prob = [], []
            for fold in kfold:
                ky_pred.append(fold.pred)
                if fold.prob is not None:
                    ky_prob.append(fold.prob)

                models.append(fold.estimator)
                scalers.append(fold.scaler)

            ky_pred = concat(ky_pred).loc[self.__y.index]
            y_pred.append(ky_pred)

            if ky_prob:
                ky_prob = concat(ky_prob).loc[self.__y.index]
                y_prob.append(ky_prob)

            for s, (f, p) in self.__scorers.items():
                fold_scorers[s].append(f(self.__y, (ky_prob if p else ky_pred)))

        y_pred = concat(y_pred, axis=1)
        y_prob = concat(y_prob, axis=1, keys=range(len(y_prob))) if y_prob else None

        res = {}
        for s, _v in fold_scorers.items():
            m, v = mean(_v), sqrt(var(_v))
            c = (1 if s == 'rmse' else -1) * m + self.__disp_coef * v
            res.update({s: m, '$%s' % s: c, '%s_var' % s: v})

        return FitContainer(models=models, scalers=scalers, params=fitparams, pred=y_pred, prob=y_prob, scores=res)

    def __cv_domain_fit(self, fitparams):
        """ NEED IMPLEMENT!
        parallel = Parallel(n_jobs=self.__n_jobs)
        folds = parallel(delayed(_kfold)(self.__domain_class, fitparams, self.__x.iloc[train], self.__y.iloc[train],
                                         self.__x.iloc[test], self.__y.iloc[test], self.__normalize)
                         for i in range(repetitions) for train, test in self.__cv_naive(i))
        #  street magic. split folds to repetitions
        for kfold in zip(*[iter(folds)] * self.__nfold):
            domains.append(kfold)
        """
        raise Exception('NOT IMPLEMENTED')

    def __domain_fit(self):
        models, scalers, y_pred, y_prob = [], [], [], []
        parallel = Parallel(n_jobs=self.__n_jobs)
        folds = parallel(delayed(_kfold)(self.__domain_class, {}, self.__x.iloc[train], self.__y.iloc[train],
                                         self.__x.iloc[test], self.__domain_normalize) for train, test in self.__cv())

        for kfold in zip(*[iter(folds)] * self.__nfold):
            ky_pred, ky_prob = [], []
            for fold in kfold:
                ky_pred.append(fold.pred)
                if fold.prob is not None:
                    ky_prob.append(fold.prob)

                models.append(fold.estimator)
                scalers.append(fold.scaler)

            ky_pred = concat(ky_pred).loc[self.__y.index]
            y_pred.append(ky_pred)

            if ky_prob:
                ky_prob = concat(ky_prob).loc[self.__y.index]
                y_prob.append(ky_prob)

        y_pred = concat(y_pred, axis=1)
        y_prob = concat(y_prob, axis=1, keys=range(len(y_prob))) if y_prob else None

        return FitContainer(models=models, scalers=scalers, params=None, pred=y_pred, prob=y_prob, scores=None)

    def __cv_naive(self):
        """ shuffling method for CV. may be smartest.
        """
        setindexes = arange(len(self.__y.index))
        for i in range(self.__repetitions):
            shuffled = shuffle(setindexes, random_state=i)
            for train, test in KFold(n_splits=self.__nfold).split(setindexes):
                yield shuffled[train], shuffled[test]

    @staticmethod
    def __dive(fit, param, md, di):
        tmp = {}
        for i, j in fit.items():
            if di[i] < md and not param[i][j]:
                tmp[i] = param[i]
            else:
                tmp[i] = param[i][j]
        if all(x for x in tmp.values()):
            return tmp

    @classmethod
    def __inception(cls, fit_params):
        max_params, fitparams, depindex, maxdep = 0, [], [], []
        for param in fit_params:
            md, mp, di, pr = 0, 1, {}, {}
            for i, j in param.items():
                mp *= len(j)
                pr[i], di[i] = cls.__split_range(j)
                if di[i] > md:
                    md = di[i]
            depindex.append(di)
            maxdep.append(md)
            fitparams.append(pr)
            max_params += mp
        return depindex, maxdep, fitparams, max_params

    @classmethod
    def __split_range(cls, param, dep=0):
        tmp = {}
        fdep = dep
        stepindex = list(range(0, len(param), round(len(param)/10) or 1))
        stepindex.insert(0, -1)
        stepindex.append(len(param))
        for i, j, k in zip(stepindex, stepindex[1:], stepindex[2:]):
            tmp[param[j]], tmpd = cls.__split_range(param[i + 1:j] + param[j + 1:k], dep=dep + 1)
            if tmpd > fdep:
                fdep = tmpd
        return tmp, fdep

    __scorers_dict = dict(rmse=(rmse, False), r2=(r2_score, False), kappa=(cohen_kappa_score, False),
                          acc=(accuracy, False))
