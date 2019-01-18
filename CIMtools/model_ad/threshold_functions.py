# -*- coding: utf-8 -*-
#
#
#  Copyright 2019 Assima Rakhimbekova <asima.astana@outlook.com>
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
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, see <https://www.gnu.org/licenses/>.
#
import numpy as np
from sklearn.utils import safe_indexing
from CIMtools.model_ad.reaction_type_control import ReactionTypeControl
from sklearn.metrics import mean_squared_error, confusion_matrix
from sklearn.model_selection import KFold, RepeatedKFold
from CIMtools.model_ad.leverage import Leverage


seed = 1
kf = KFold(n_splits=5, random_state=seed, shuffle=True)
rkf = RepeatedKFold(n_splits=5, random_state=seed, n_repeats=10)


def outputY_AD(ad, X, y, metric, reg_model, ad_model=None, y_clf=None, data=None):
    """
    Нужен для оптимизации гиперпараметров АД
    Тип как мой gridSerachCV
    Ниже есть налоагичная функция outputY_AD__ отличие в том что в outputY_AD функции в качестве у подаются истинные
    значения у (которые являются external training set) и тут строится модель на internal training set, предсказывается
    internal test set
    А в функции outputY_AD__ в качестве у подается 2d array (где в первом столбце предсказанные значения external training set
    а в столбце два - true y for external training set

    Пример
    --------
    path = '/home/asima/asima/FR/datasets' # путь до датасетов
    name = 'SN2_article_2018_final++.rdf' # название датасета 1 - обучающая выборка  Нужно будет еще для названия пиклов
    constant_name = 'logK' # Таутомеры - 'tabulated_constant'
    file = '{}/{}'.format(path, name)
    environ["PATH"]+=":/home/asima" # путь до фрагментора
    for score in ['ba', 'rmse']:
        print(score)
        AD_gpr = [[] for _ in range(3)] # для gpr нужен свой список так как предсказанные значения Y у GPR свои
        if score == 'ba':
            AD_res = [[] for _ in range(36)] # для каждого AD будет сохраняться значения AD в виде True/False. чтобы
            # потом посчитать общую таблицу
        else:
            AD_res = [[] for _ in range(40)] # 4  метода считаются только при rmse (в самом конце)
        Y = []  # в этом листе сохраним все значения у для расчета полной таблицы
        Y_gpr = []  # тоже самое только для гаусовых процессов так как это отдельная регрисионная модель
        seed = 1  # for random state
        kf = KFold(n_splits=5, random_state=seed, shuffle=True)  # for GridSearchCV
        rkf = RepeatedKFold(n_splits=5, random_state=seed, n_repeats=10)  # for rkf_method (it is one type of AD)
        reactions = RDFread(file).read()

        for train_index_ext, test_index_ext in kf.split(reactions):  # external set
            reactions_train = safe_indexing(reactions, train_index_ext)
            reactions_test = safe_indexing(reactions, test_index_ext)
            data = initial_data(reactions_train, constant_name, reactions_test, constant_name)
            X_train = data[0][0] # np.array
            Y_train = data[0][1] # list
            y_train_gpr = data[0][2][:, 0] # np.array
            X_test = data[1][0]  # np.array
            Y_test = data[1][1]  # list
            y_test_gpr = data[1][2][:, 0]  # np.array
            # ------------------------------------- Строим регресионные модели ---------------------------------------------
            est = GridSearchCV(RandomForestRegressor(random_state=seed),
                               {'max_features': [0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 'auto', 'log2', None], 'n_estimators': [500]},
                               cv=kf, verbose=1, scoring='neg_mean_squared_error', n_jobs=4).fit(X_train, Y_train)  # модель
            Y_predicted = cross_val_predict(est.best_estimator_, X_train, Y_train, cv=kf, verbose=1)
            #---------вот тут получаем Y_pr_ts который и будет входным данным для функции outputY_AD__ ----------------
            Y_pr_ts = np.column_stack((Y_predicted, Y_train))
            ....

    Поэтому поидеи можно сократить объем вычисления подавая Y_pr_ts в качестве external training set и пользуясь
    не этой функцией а следующей. Но так как надо понять логику того как находятся гиперпараметры ад оставила эту функцию

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
            y_pred = abs(y_test - y_pred) <= 3 * np.sqrt(mean_squared_error(y_test, y_pred))  # для расчета ba
        Y_pred.extend(y_pred)
        Y_true.extend(y_test)

        if ad == 'rtc':
            ad_model.fit(data, index=train_index)
            AD.extend(ad_model.predict(data, index=test_index))
        elif ad == 'leverage':
            ad_model = Leverage().fit(x_train)
            AD.extend(ad_model.predict_proba(x_test))
        elif ad == 'simdist':
            AD.extend(data_test)
        elif ad == 'rfr_var':
            ad_model.fit(x_train, y_train)
            AD.extend(ad_model.predict_proba(x_test))
        elif ad == '2cc':
            ad_model.fit(x_train, y_train_clf)
            AD.extend(ad_model.predict_proba(x_test)[:, 0])
        elif ad == 'gpr':
            AD.extend(y_var)
    return np.array(Y_true), np.array(Y_pred), np.array(AD)


def my_helper_1(metric, y_test):
    """
    эта функция нужна если пользоваться функцией outputY_AD__

    :param metrics: по чему оптимизируем гиперпараметры ад
    :param y_test: Y_pred_test: np.array 2D, shape: (n,2)
     First column Y_predicted by RandomForest Regressor for internal set(x_test)
     Second column Y_test, experimental values y_test
     Y_pred = we need to obtain np.array consist of only True and False, so we compare first and second
     column and compare this difference with 2RMSE or 3RMSE
    :return: y_pred
    """
    if metric == 'r2':
        y_pred = y_test[:, 0]
    elif metric == 'ba':
        y_pred = abs(y_test[:, 0] - y_test[:, 1]) <= 3 * np.sqrt(mean_squared_error(y_test[:, 1], y_test[:, 0]))
    elif metric == 'auc':
        y_pred = -abs(y_test[:, 0] - y_test[:, 1])
    elif metric == 'rmse':
        y_pred = y_test[:, 0]
    elif metric == 'coverage':
        y_pred = y_test[:, 0]
    return y_pred


def outputY_AD__(ad, X, y, metric, reg_model, ad_model=None, y_clf=None, data=None):
    """
    Нужен для оптимизации гиперпараметров АД
    Тип как мой gridSerachCV
    Ниже есть налоагичная функция outputY_AD__ отличие в том что в outputY_AD функции в качестве у подаются истинные
    значения у (которые являются external training set) и тут строится модель на internal training set, предсказывается
    internal test set
    А в функции outputY_AD__ в качестве у подается 2d array (где в первом столбце предсказанные значения external training set
    а в столбце два - true y for external training set

    Пример
    --------
    path = '/home/asima/asima/FR/datasets' # путь до датасетов
    name = 'SN2_article_2018_final++.rdf' # название датасета 1 - обучающая выборка  Нужно будет еще для названия пиклов
    constant_name = 'logK' # Таутомеры - 'tabulated_constant'
    file = '{}/{}'.format(path, name)
    environ["PATH"]+=":/home/asima" # путь до фрагментора
    for score in ['ba', 'rmse']:
        print(score)
        AD_gpr = [[] for _ in range(3)] # для gpr нужен свой список так как предсказанные значения Y у GPR свои
        if score == 'ba':
            AD_res = [[] for _ in range(36)] # для каждого AD будет сохраняться значения AD в виде True/False. чтобы
            # потом посчитать общую таблицу
        else:
            AD_res = [[] for _ in range(40)] # 4  метода считаются только при rmse (в самом конце)
        Y = []  # в этом листе сохраним все значения у для расчета полной таблицы
        Y_gpr = []  # тоже самое только для гаусовых процессов так как это отдельная регрисионная модель
        seed = 1  # for random state
        kf = KFold(n_splits=5, random_state=seed, shuffle=True)  # for GridSearchCV
        rkf = RepeatedKFold(n_splits=5, random_state=seed, n_repeats=10)  # for rkf_method (it is one type of AD)
        reactions = RDFread(file).read()

        for train_index_ext, test_index_ext in kf.split(reactions):  # external set
            reactions_train = safe_indexing(reactions, train_index_ext)
            reactions_test = safe_indexing(reactions, test_index_ext)
            data = initial_data(reactions_train, constant_name, reactions_test, constant_name)
            X_train = data[0][0] # np.array
            Y_train = data[0][1] # list
            y_train_gpr = data[0][2][:, 0] # np.array
            X_test = data[1][0]  # np.array
            Y_test = data[1][1]  # list
            y_test_gpr = data[1][2][:, 0]  # np.array
            # ------------------------------------- Строим регресионные модели ---------------------------------------------
            est = GridSearchCV(RandomForestRegressor(random_state=seed),
                               {'max_features': [0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 'auto', 'log2', None], 'n_estimators': [500]},
                               cv=kf, verbose=1, scoring='neg_mean_squared_error', n_jobs=4).fit(X_train, Y_train)  # модель
            Y_predicted = cross_val_predict(est.best_estimator_, X_train, Y_train, cv=kf, verbose=1)
            #---------вот тут получаем Y_pr_ts который и будет входным данным для функции outputY_AD__ ----------------
            Y_pr_ts = np.column_stack((Y_predicted, Y_train))
            ....

    Поэтому поидеи можно сократить объем вычисления подавая Y_pr_ts в качестве external training set и пользуясь
    не этой функцией а следующей. Но так как надо понять логику того как находятся гиперпараметры ад оставила эту функцию

    :param ad: название области применимости модели: rfr_var, rtc, leverage, gpr, simdist, 2cc
    :param X: X_train (обучающая выборка)
    :param y: Y_pred_test: np.array 2D, shape: (n,2)
     First column Y_predicted by RandomForest Regressor for internal set(x_test)
     Second column Y_test, experimental values y_test
     Y_pred = we need to obtain np.array consist of only True and False, so we compare first and second
     column and compare this difference with 2RMSE or 3RMSE
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
            y_train_clf = safe_indexing(y_clf, train_index)
            y_test_clf = safe_indexing(y_clf, test_index)
        if data is not None:
            data_train = safe_indexing(data, train_index)
            data_test = safe_indexing(data, test_index)

        if ad == 'gpr':
            reg_model.fit(x_train, y_train)
            y_pred, y_var = reg_model.predict(x_test, return_std=True)
            if metric == 'ba':
                y_pred = abs(y_test - y_pred) <= 3 * np.sqrt(mean_squared_error(y_test, y_pred))
                Y_pred.extend(y_pred)
            else:
                Y_pred.extend(y_pred)
            Y_true.extend(y_test)
        else:
            y_pred = my_helper_1(metric, y_test)
            Y_pred.extend(y_pred)
            Y_true.extend(y_test[:, 1])

        if ad == 'rtc':
            ad_model.fit(data, index=train_index)
            AD.extend(ad_model.predict(data, index=test_index)) # list
        elif ad == 'leverage':
            ad_model = Leverage().fit(x_train)
            AD.extend(ad_model.predict_proba(x_test)) # array
        elif ad == 'simdist':
            AD.extend(data_test) # array
        elif ad == 'rfr_var':
            ad_model.fit(x_train, y_train[:, 1])
            AD.extend(ad_model.predict_proba(x_test)) # array to list
        elif ad == '2cc':
            ad_model.fit(x_train, y_train_clf)
            AD.extend(ad_model.predict_proba(x_test)[:, 0]) # array
        elif ad == 'gpr':
            AD.extend(y_var) # array
    return np.array(Y_true), np.array(Y_pred), np.array(AD)


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
            RMSE_AD = np.sqrt((sum(map(lambda x: (((x[0] - x[1]) ** 2) * x[2]), zip(Y_pred, Y_true, AD_pred_n)))) / s_n)
        else:
            RMSE_AD = 0

        if s_out_n:
            RMSE_AD_out_n = np.sqrt(
                (sum(map(lambda x: (((x[0] - x[1]) ** 2) * x[2]), zip(Y_pred, Y_true, AD_pred_out_n)))) / s_out_n)
        else:
            RMSE_AD_out_n = 0
        val = RMSE_AD_out_n - RMSE_AD
    return val


def threshold(ad, X, y, metric, reg_model, ad_model=None, envs=None, y_clf=None, data=None):
    """
    Функция для подбора оптимаьных гиперпараметров

    :param ad: название ад (в виде str like 'rtc', 'gpr, 'leverage', 'simdist')
    :param X: x for external training set
    :param y: Y_pred_test: np.array 2D, shape: (n,2)
     First column Y_predicted by RandomForest Regressor for internal set(x_test)
     Second column Y_test, experimental values y_test
     Y_pred = we need to obtain np.array consist of only True and False, so we compare first and second
     column and compare this difference with 2RMSE or 3RMSE
     or
    just y but NB! тогда надо исправить код! и пользоватья функцией outputY_AD
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
            ad_model = ReactionTypeControl(env=env, hyb=True)
            results = outputY_AD__(ad=ad, X=X, y=y, metric=metric, reg_model=reg_model, ad_model=ad_model, y_clf=None,
                                   data=data)
            val = function(metric=metric, Y_true=results[0], Y_pred=results[1], AD=results[2])
            if val >= min_hparam['score']:
                min_hparam['score'] = val
                min_hparam['z'] = env
        return min_hparam
    elif ad == 'leverage':
        results = outputY_AD__(ad=ad, X=X, y=y, metric=metric, reg_model=reg_model)
    elif ad == 'simdist':
        results = outputY_AD__(ad=ad, X=X, y=y, metric=metric, reg_model=reg_model, ad_model=None, y_clf=None, data=data)
    else:
        results = outputY_AD__(ad=ad, X=X, y=y, metric=metric, reg_model=reg_model, ad_model=ad_model, y_clf=y_clf)
    AD_ = np.unique(results[2])
    for z in AD_:
        AD_new = results[2] <= z
        val = function(metric=metric, Y_true=results[0], Y_pred=results[1], AD=AD_new)
        if val >= min_hparam['score']:
            min_hparam['score'] = val
            min_hparam['z'] = z
    return min_hparam


def my_score_ba(Y_pred_test, AD_pred):
    """
    Function for gridsearchsv, scoring function
    it must put for make_scorer() function
     :param Y_pred_test: np.array 2D, shape: (n,2)
     First column Y_predicted by RandomForest Regressor for internal set(x_test)
     Second column Y_test, experimental values y_test
     Y_pred = we need to obtain np.array consist of only True and False, so we compare first and second
     column and compare this difference with 2RMSE or 3RMSE
     :param AD_pred: results for 1-SVM (+1 or -1), IF(+1,0) as np.array with shape (n,)
     Then we need to perform AD_pred consists only True and False
     :return: ba
     """
    Y_pred = abs(Y_pred_test[:, 0] - Y_pred_test[:, 1]) <= 3 * np.sqrt(mean_squared_error(Y_pred_test[:, 1],
                                                                                          Y_pred_test[:, 0]))
    AD_pred_n = (AD_pred == 1)
    tn, fp, fn, tp = confusion_matrix(Y_pred, AD_pred_n).ravel()
    sen_TPR = (float(tp)) / (tp + fn)
    spc_TNR = (float(tn) / (tn + fp))
    ba = ((sen_TPR + spc_TNR) / 2)
    return ba


def my_score_rmse(Y_pred_test, AD_pred):
    """
    Function for gridsearchsv, scoring function
    it must put for make_scorer() function
     :param Y_pred_test: np.array 2D, shape: (n,2)
     First column Y_predicted by RandomForest Regressor for internal set(x_test)
     Second column Y_test, experimental values y_test
     Y_pred = we need to obtain np.array consist of only True and False, so we compare first and second
     column and compare this difference with 2RMSE or 3RMSE
     :param AD_pred: results for 1-SVM (+1 or -1), IF(+1,0) as np.array with shape (n,)
     Then we need to perform AD_pred consists only True and False
     :return: diff_RMSE
     """
    AD_pred_n = (AD_pred == 1)
    AD_pred_out = AD_pred_n - 1
    AD_pred_out_n = (AD_pred_out == -1)
    s_n = sum(AD_pred_n)
    s_out_n = sum(AD_pred_out_n)
    if s_n:
        RMSE_AD = np.sqrt((sum(map(lambda x: (((x[0] - x[1]) ** 2) * x[2]), zip(Y_pred_test[:, 0], Y_pred_test[:, 1],
                                                                             AD_pred_n)))) / s_n)
    else:
        RMSE_AD = 0

    if s_out_n:
        RMSE_AD_out_n = np.sqrt((sum(map(lambda x: (((x[0] - x[1]) ** 2) * x[2]), zip(Y_pred_test[:, 0], Y_pred_test[:, 1],
                                                                                   AD_pred_out_n)))) / s_out_n)
    else:
        RMSE_AD_out_n = 0
    diff_RMSE = RMSE_AD_out_n - RMSE_AD
    return diff_RMSE


def UniqueZ(Y_ext, Y_pred, AD, metric, more):
    """
    Эта функция используется только в consensus_pred_var. Также для поиска оптимальных гиперпараметров
    :param Y_ext:
    :param Y_pred:
    :param AD:
    :param metric:
    :param more:
    :return:
    """
    min_hparam = {'z': 0, 'score': 0}
    AD_ = np.unique(AD)
    for z in AD_:
        if more == '>':
            AD_new = AD > z
        elif more == '>=':
            AD_new = AD >= z
        elif more == '<':
            AD_new = AD < z
        elif more == '<=':
            AD_new = AD <= z
        val = function(metric=metric, Y_true=Y_ext, Y_pred=Y_pred, AD=AD_new)
        if val >= min_hparam['score']:
            min_hparam['score'] = val
            min_hparam['z'] = z
    return min_hparam
