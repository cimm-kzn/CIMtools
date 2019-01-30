# -*- coding: utf-8 -*-
#
#  Copyright 2019 Assima Rakhimbekova <asima.astana@outlook.com>
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
from CIMtools.applicability_domain.reaction_type_control import ReactionTypeControl
from CIMtools.applicability_domain.leverage import Leverage
from numpy import array, unique, sqrt
from sklearn.utils import safe_indexing
from sklearn.metrics import mean_squared_error, confusion_matrix
from sklearn.model_selection import KFold, RepeatedKFold


kf = KFold(n_splits=5, random_state=1, shuffle=True)
rkf = RepeatedKFold(n_splits=5, random_state=1, n_repeats=10)

def outputY_AD(ad, X, y, metric, reg_model, ad_model=None, y_clf=None, data=None):
    """
    Нужен для оптимизации гиперпараметров АД
    Тип как мой gridSerachCV
    :param ad: название области применимости модели: rfr_var, rtc, leverage, gpr, simdist, 2cc
    :param X: X_train (обучающая выборка)
    :param y: y_train (обучающая выборка)
    :param metric: метрика, по которой происходит оптимизация: 'ba', 'rmse', 'r2_score', 'auc', 'coverage' и тд и тп
    :param reg_model: регресионная модель, для построения модели и предсказания свойств (y)
    :param ad_model: модель АД для предсказания вероятности того на сколько объект принадлежит АД (все кроме)
     либо просто True/False - 'rtc'
    :param y_clf: у для RandomForestClassifier в виде (True/False) True - те ошибка предсказания которых меньше 3 RMSE
    :param data: нужен для других случаев. Например, Для simdist - это расстояния до соседа и при этом это и есть АД.
    Для RTC - это rdf file
    :return: np.array(Y_true), np.array(Y_pred), np.array(AD)
    """
    Y_pred = []
    Y_true = []
    AD = []
    for train_index, test_index in kf.split(X):
        x_train = safe_indexing(X, train_index)
        x_test = safe_indexing(X, test_index)
        y_train = safe_indexing(y, train_index)
        y_test = safe_indexing(y, test_index)
        if y_clf is not None:
            y_train_clf = safe_indexing(y_clf, train_index) # for RandomForestClassifier
            y_test_clf = safe_indexing(y_clf, test_index)
        if data is not None:
            data_train = safe_indexing(data, train_index) # for simdist, rtc
            data_test = safe_indexing(data, test_index)

        reg_model.fit(x_train, y_train) # обучаем модель будь то RandomForestRegressor or GPR
        if ad == 'gpr':
            y_pred, y_var = reg_model.predict(x_test, return_std=True) # y_var - это дисперсия предсказанных значений
        else:
            y_pred = reg_model.predict(x_test) # просто предсказанные значения

        if metric == 'ba':
            y_pred = abs(y_test - y_pred) <= 3 * sqrt(mean_squared_error(y_test, y_pred))  # для расчета ba
        Y_pred.extend(y_pred)
        Y_true.extend(y_test)

        if ad == 'rtc':
            ad_model.fit(data_train)
            AD.extend(ad_model.predict(data_test))
        elif ad == 'leverage':
            ad_model = Leverage().fit(x_train)
            AD.extend(ad_model.predict_proba(x_test))
    return array(Y_true), array(Y_pred), array(AD)

def function(metric, Y_true, Y_pred, AD):
    """
     для расчета скоринговой функции
    :param metric: по чему оптимизируем
    :param Y_true: истиные значения у for external training set 1d array)
    :param Y_pred: предсказанные значения у for external training set (1d array)
    :param AD: 1d array)
    :return: ba or diff_RMSE (is called val)
    """
    if metric == 'ba':
        AD_pred_n = (AD == 1)
        tn, fp, fn, tp = confusion_matrix(Y_pred, AD_pred_n).ravel()
        sen_TPR = (float(tp))/(tp+fn)
        spc_TNR = (float(tn)/(tn+fp))
        val = ((sen_TPR+spc_TNR)/2)
    else:
        AD_pred_n = (AD == 1)
        AD_pred_out = AD_pred_n - 1
        AD_pred_out_n = (AD_pred_out == -1)

        s_n = sum(AD_pred_n)
        s_out_n = sum(AD_pred_out_n)
        if s_n:
            RMSE_AD = sqrt((sum(map(lambda x: (((x[0] - x[1]) ** 2) * x[2]), zip(Y_pred, Y_true, AD_pred_n)))) / s_n)
        else:
            RMSE_AD = 0

        if s_out_n:
            RMSE_AD_out_n = sqrt(
                (sum(map(lambda x: (((x[0] - x[1]) ** 2) * x[2]), zip(Y_pred, Y_true, AD_pred_out_n)))) / s_out_n)
        else:
            RMSE_AD_out_n = 0
        val = RMSE_AD_out_n - RMSE_AD
    return val


def threshold(ad, X, y, metric, reg_model, ad_model=None, envs=None, y_clf=None, data=None):
    """
    Функция для подбора оптимальных гиперпараметров

    :param ad: название ад (в виде str like 'rtc', 'gpr, 'leverage', 'simdist')
    :param X: x for external training set
    :param y: y
    :param metric: по чему оптимизировать ('ba', 'rmse')
    :param reg_model: регресионная модель (в моем случае RandomForestRegression, которая была ранее с оптимизированна)
    :param ad_model: для 'rfr_var'
    :param envs: for 'rtc'
    :param y_clf: for '2cc
    :param data: for 'rtc'
    :return: гиперпараметр (причем у меня пока написано для таких методов, у которых оптимизируется только один гиперпараметр
    да так чтобы он максимизировал скоринговую функцию
    """
    min_hparam = {'z': 0, 'score': 0}
    if ad == 'rtc':
        for env in envs:
            ad_model = ReactionTypeControl(env=env)
            results = outputY_AD(ad=ad, X=X, y=y, metric=metric, reg_model=reg_model, ad_model=ad_model, data=data)
            val = function(metric=metric, Y_true=results[0], Y_pred=results[1], AD=results[2])
            if val >= min_hparam['score']:
                min_hparam['score'] = val
                min_hparam['z'] = env
        return min_hparam
    elif ad == 'leverage':
        results = outputY_AD(ad=ad, X=X, y=y, metric=metric, reg_model=reg_model)
    elif ad == 'simdist':
        results = outputY_AD(ad=ad, X=X, y=y, metric=metric, reg_model=reg_model, data=data)
    else:
        results = outputY_AD(ad=ad, X=X, y=y, metric=metric, reg_model=reg_model, ad_model=ad_model)
    AD_ = unique(results[2])
    for z in AD_:
        AD_new = results[2] <= z
        val = function(metric=metric, Y_true=results[0], Y_pred=results[1], AD=AD_new)
        if val >= min_hparam['score']:
            min_hparam['score'] = val
            min_hparam['z'] = z
    return min_hparam
