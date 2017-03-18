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
from collections import defaultdict
from copy import deepcopy
from itertools import product
from math import sqrt, ceil
from numpy import inf, mean, var, arange
from operator import lt, le
from pandas import DataFrame, Series, concat
from shutil import rmtree
from sklearn.cross_validation import KFold
from sklearn.externals.joblib import Parallel, delayed
from sklearn.metrics import mean_squared_error, r2_score, cohen_kappa_score, accuracy_score
from sklearn.preprocessing import MinMaxScaler
from sklearn.utils import shuffle
from tempfile import mkdtemp
from typing import Iterable


class Score(dict):
    def __init__(self, *args, **kwargs):
        dict.__init__(self, *args, **kwargs)

    def __comparer(self, other, op):
        return all(op(y, other[x]) for x, y in self.items())

    def __lt__(self, other):
        return self.__comparer(other, lt)

    def __le__(self, other):
        return self.__comparer(other, le)


def _kfold(est, x, y, train, test, svmparams, normalize, box):
    x_train, y_train = x.iloc[train], y.iloc[train]
    x_test, y_test = x.iloc[test], y.iloc[test]
    x_min = x_train.min().loc[box]
    x_max = x_train.max().loc[box]
    y_min = y_train.min()
    y_max = y_train.max()

    x_ad = ((x_test.loc[:, box] - x_min).min(axis=1) >= 0) & ((x_max - x_test.loc[:, box]).min(axis=1) >= 0)

    if normalize:
        normal = MinMaxScaler()
        x_train = DataFrame(normal.fit_transform(x_train), columns=x_train.columns)
        x_test = DataFrame(normal.transform(x_test), columns=x_train.columns)
    else:
        normal = None

    model = est(**svmparams)
    model.fit(x_train, y_train)
    y_pred = Series(model.predict(x_test), index=y_test.index)

    y_ad = (y_pred >= y_min) & (y_pred <= y_max)

    output = dict(model=model, normal=normal, x_min=x_min, x_max=x_max, y_min=y_min, y_max=y_max,
                  y_test=y_test, y_pred=y_pred, x_ad=x_ad, y_ad=y_ad)

    if hasattr(model, 'predict_proba'):
        output['y_prob'] = DataFrame(model.predict_proba(x_test), index=y_test.index, columns=model.classes_)

    return output


def _rmse(y_test, y_pred):
    return sqrt(mean_squared_error(y_test, y_pred))


def _accuracy(y_test, y_pred):
    return accuracy_score(y_test, y_pred, normalize=True)


def _iap(y_test, y_prob):
    res = defaultdict(list)
    for sid, s_prob in y_prob.groupby(level='structure'):
        s_test = y_test.xs(sid, level='structure', drop_level=False)
        for col in s_prob.columns:
            in_class = s_test.loc[s_test == col].index
            out_class = s_test.loc[s_test != col].index

            in_class_dominance = sum(s_prob.loc[i][col] > s_prob.loc[o][col] for i, o in product(in_class, out_class))
            possible_combo = len(in_class) * len(out_class)
            if possible_combo:
                res[col].append(in_class_dominance / possible_combo)
    res = Score({x: sum(y) / len(y) for x, y in res.items()})
    return res


class BaseModel(ABC):
    def __init__(self, descriptors_generator, structures, workpath='.', nfold=5, repetitions=1, rep_boost=100,
                 dispcoef=0, fit='rmse', scorers=('rmse', 'r2'), normalize=False, n_jobs=2, **kwargs):

        _scorers = dict(rmse=(_rmse, False), r2=(r2_score, False),
                        kappa=(cohen_kappa_score, False), acc=(_accuracy, False), iap=(_iap, True))
        self.__model = {}

        self.__generator = descriptors_generator
        self.set_work_path(workpath)

        self.__nfold = nfold
        self.__repetitions = repetitions
        self.__rep_boost = ceil(repetitions * (rep_boost % 100) / 100) or repetitions
        self.__n_jobs = n_jobs

        self.__normalize = normalize
        self.__disp_coef = dispcoef
        self.__scorers = {x: _scorers[x] for x in scorers if x in _scorers}
        self.__fit_score = 'C' + (fit if fit in scorers else scorers[0])
        self.__score_reporter = '\n'.join(['{0} +- variance = %({0})s +- %(v{0})s'.format(i) for i in self.__scorers])

        print("Descriptors generation start")
        xy = descriptors_generator.get(structures, **kwargs)
        self.__x = xy['X']
        self.__y = xy['Y']
        self.__box = xy.get('BOX', xy['X'].columns)
        print("Descriptors generated")

        self.__cross_val()
        self.delete_work_path()

    @abstractmethod
    def prepare_params(self, param):
        pass

    @property
    @abstractmethod
    def fit_params(self) -> Iterable:
        pass

    @property
    @abstractmethod
    def estimator(self):
        pass

    def set_work_path(self, workpath):
        self.__workpath = mkdtemp(dir=workpath)
        if hasattr(self.__generator, 'set_work_path'):
            self.__generator.set_work_path(self.__workpath)

    def delete_work_path(self):
        rmtree(self.__workpath)

    def get_model_stats(self):
        stat = {x: self.__model[x] for x in self.__scorers}
        stat.update({'%s_var' % x: self.__model['v%s' % x] for x in self.__scorers})

        stat.update(dict(fitparams=self.__model['params'], repetitions=self.__repetitions,
                         nfolds=self.__nfold, normalize=self.__normalize, dragostolerance=sqrt(self.__y.var())))
        return stat

    def get_fit_predictions(self):
        output = dict(property=self.__y, prediction=self.__model['y_pred'], y_domain=self.__model['y_ad'],
                      domain=self.__model['x_ad'])
        if 'y_prob' in self.__model:
            output['probability'] = self.__model['y_prob']
        return output

    def get_models(self):
        return self.__model['models']

    def get_descriptors_generator(self):
        return self.__generator

    def __split_range(self, param, dep=0):
        tmp = {}
        fdep = dep
        stepindex = list(range(0, len(param), round(len(param)/10) or 1))
        stepindex.insert(0, -1)
        stepindex.append(len(param))
        for i, j, k in zip(stepindex, stepindex[1:], stepindex[2:]):
            tmp[param[j]], tmpd = self.__split_range(param[i + 1:j] + param[j + 1:k], dep=dep + 1)
            if tmpd > fdep:
                fdep = tmpd
        return tmp, fdep

    def __cross_val(self):
        fitparams = deepcopy(self.fit_params)
        fcount = 0
        depindex = []
        maxdep = []
        print('list of fit params:')
        print(DataFrame(list(fitparams)))
        for param in fitparams:
            di = {}
            md = 0
            for i in param:
                if i != 'kernel':
                    param[i], di[i] = self.__split_range(param[i])
                    if di[i] > md:
                        md = di[i]
            depindex.append(di)
            maxdep.append(md)

        print('========================================\n'
              'Y mean +- variance = %s +- %s\n'
              '  max = %s, min = %s\n'
              '========================================' %
              (self.__y.mean(), sqrt(self.__y.var()), self.__y.max(), self.__y.min()))

        bestmodel = badmodel = dict(model=None, Cr2=inf, Crmse=inf, Ckappa=inf, Cba=inf, Ciap=inf)
        for param, md, di in zip(fitparams, maxdep, depindex):
            var_kern_model = badmodel
            while True:
                var_param_model = badmodel
                tmp = self.prepare_params(param)
                for i in tmp:
                    fcount += 1
                    print('%d: fit model with params:' % fcount, i)
                    fittedmodel = self.__fit(i, self.__rep_boost)
                    print(self.__score_reporter % fittedmodel)
                    if fittedmodel[self.__fit_score] < var_param_model[self.__fit_score]:
                        var_param_model = fittedmodel

                if var_param_model[self.__fit_score] < var_kern_model[self.__fit_score]:
                    var_kern_model = var_param_model
                    tmp = {}
                    for i, j in var_kern_model['params'].items():
                        if i in ('kernel', 'probability'):
                            tmp[i] = j
                        elif di[i] < md and not param[i][j]:
                            tmp[i] = param[i]
                        else:
                            tmp[i] = param[i][j]
                    param = tmp
                else:
                    break
            if var_kern_model[self.__fit_score] < bestmodel[self.__fit_score]:
                bestmodel = var_kern_model

        if self.__repetitions > self.__rep_boost:
            bestmodel = self.__fit(bestmodel['params'], self.__repetitions)
        print('========================================\n' +
              ('SVM params %(params)s\n' + self.__score_reporter) % bestmodel)
        print('========================================\n%s variants checked' % fcount)
        self.__model = bestmodel

    def __fit(self, fitparams, repetitions):
        models, y_pred, y_prob, y_ad, x_ad = [], [], [], [], []
        fold_scorers = defaultdict(list)
        parallel = Parallel(n_jobs=self.__n_jobs)
        kf = list(KFold(len(self.__y), n_folds=self.__nfold))
        setindexes = arange(len(self.__y.index))
        folds = parallel(delayed(_kfold)(self.estimator, self.__x, self.__y, s[train], s[test],
                                         fitparams, self.__normalize, self.__box)
                         for s in (self.__shuffle(setindexes, i) for i in range(repetitions))
                         for train, test in kf)

        #  street magic. split folds to repetitions
        for kfold in zip(*[iter(folds)] * self.__nfold):
            ky_pred, ky_prob, ky_ad, kx_ad = [], [], [], []
            for fold in kfold:
                ky_pred.append(fold.pop('y_pred'))
                ky_ad.append(fold.pop('y_ad'))
                kx_ad.append(fold.pop('x_ad'))

                if 'y_prob' in fold:
                    ky_prob.append(fold.pop('y_prob'))

                fold.pop('y_test')
                models.append(fold)

            ky_pred = concat(ky_pred).loc[self.__y.index]
            ky_ad = concat(ky_ad).loc[self.__y.index]
            kx_ad = concat(kx_ad).loc[self.__y.index]

            if ky_prob:
                ky_prob = concat(ky_prob).loc[self.__y.index].fillna(0)
                y_prob.append(ky_prob)

            for s, (f, p) in self.__scorers.items():
                fold_scorers[s].append(f(self.__y, (ky_prob if p else ky_pred)))

            y_pred.append(ky_pred)
            y_ad.append(ky_ad)
            x_ad.append(kx_ad)

        y_pred = concat(y_pred, axis=1)
        y_ad = concat(y_ad, axis=1)
        x_ad = concat(x_ad, axis=1)

        res = dict(models=models, params=fitparams, y_pred=y_pred, y_ad=y_ad, x_ad=x_ad)
        if y_prob:
            res['y_prob'] = concat(y_prob, axis=1, keys=range(len(y_prob)))

        for s, _v in fold_scorers.items():
            if isinstance(_v[0], Score):
                m, v, c = Score(), Score(), Score()
                tmp = defaultdict(list)
                for _s in _v:
                    for k, val in _s.items():
                        tmp[k].append(val)
                for k, val in tmp.items():
                    m[k] = mean(val)
                    v[k] = sqrt(var(val))
                    c[k] = -m[k] + self.__disp_coef * v[k]
            else:
                m, v = mean(_v), sqrt(var(_v))
                c = (1 if s == 'rmse' else -1) * m + self.__disp_coef * v
            res.update({s: m, 'C%s' % s: c, 'v%s' % s: v})

        return res

    @staticmethod
    def __shuffle(setindexes, seed):
        """ shuffling method for CV. may be smartest.
        """
        shuffled = shuffle(setindexes, random_state=seed)
        return shuffled

    def predict(self, structures, **kwargs):
        res = self.__generator.get(structures, **kwargs)
        d_x, d_ad, d_s = res['X'], res['AD'], res.get('structures')

        pred, prob, x_ad, y_ad = [], [], [], []
        for i, model in enumerate(self.__model['models']):
            x_t = DataFrame(model['normal'].transform(d_x), columns=d_x.columns) if model['normal'] else d_x

            y_p = Series(model['model'].predict(x_t), index=d_x.index)
            pred.append(y_p)

            if hasattr(model['model'], 'predict_proba'):
                y_pa = DataFrame(model['model'].predict_proba(x_t), index=d_x.index, columns=model['model'].classes_)
                prob.append(y_pa)

            y_ad.append((y_p >= model['y_min']) & (y_p <= model['y_max']))
            x_ad.append(((d_x.loc[:, self.__box] - model['x_min']).min(axis=1) >= 0) &
                        ((model['x_max'] - d_x.loc[:, self.__box]).min(axis=1) >= 0) & d_ad)

        out = dict(prediction=concat(pred, axis=1),
                   domain=concat(x_ad, axis=1), y_domain=concat(y_ad, axis=1))

        if prob:
            out['probability'] = concat(prob, axis=1, keys=range(len(prob)))

        if d_s is not None:
            out['structures'] = d_s

        return out
