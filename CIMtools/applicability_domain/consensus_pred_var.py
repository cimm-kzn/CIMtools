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
from sklearn.model_selection import RepeatedKFold, KFold
from pandas import Series, concat
from CIMtools.applicability_domain.bounding_box import Box
from CIMtools.applicability_domain.domain_selection.threshold_functions import UniqueZ
from sklearn.base import BaseEstimator
from sklearn.utils.validation import check_array, check_is_fitted
from sklearn.ensemble import RandomForestRegressor
from CIMtools.applicability_domain.reaction_type_control import ReactionTypeControl


seed = 1
kf = KFold(n_splits=5, random_state=seed, shuffle=True)
rkf = RepeatedKFold(n_splits=5, random_state=seed, n_repeats=10)


class AD_rkf(BaseEstimator):
    """
    делаем множество выборок (50 models), на каждой выборке строим модели и AD (BB+FC), причем BB вызываем в этом классе,
    то бишь внутри, а FC or RTC(deep=1) вызывааем снуражи, то есть зафитинный на external training set,
    То есть у нас список У = где для каждой реакции предсказано свойство 50 раз
    AD = 50 True/False для каждой реакции. Затем есть два пути а следовательно и два АД
    1) затем вычисляем либо средние значения предсказанный свойств и среднее значение AD - такой АД назвала consensus,
    2) либо в качестве AD дисперсия предсказанных значений У - назвала такой АД pred_var.

    Тогда получается список для каждой реакции среднее значение у, дисперсия у, среднее значение в АД. Получается в АД числа

    В качестве внутренних параметров значения отсечки. Таким образом, по согласию этих моделей определялось объект
    входит в AD или нет.

    1)словарь состоящий из индекса реакции и всех его предсказанных значений и снова словарь также с индексами
    только результатов АД
                                                    *** I ***
                                                   консенсусный
    2.1) найти средние значения для каждой строки для предсказанных значений
    2.2) процент каждой строки (каждой реакции) что она в АД
    3) Вот у нас вектор средних предсказанных значений У и его процент в АД, надо найти такую отсечку которая
    максимизирует мою метрику (R2_score or ba)
    4) Затем если вернуться к экстренал сету то каждый раз строят модель на интренал сете и для
    экстренал сета и надо также искать ВВ и вызывать FC
    5) в итоге есть вектор средних предсказанных значений и также проценты в АД и по уже найденной отсечке
    значения ниже которой будут в АД выше в не АД
    потом получаем ветор true and False и есть вектор средних предсказанных значений также рассчитыаем метрики
                                                    *** II ****
                                             дисперсия предсказаний
    В данном случае в качестве АД будет словарь также из индексов и дсперсии предскаанных значений у
    Также ищем отсечку и делаем на экстренал сете
    """
    def __init__(self, metric='rmse', reg_model=RandomForestRegressor(n_estimators=500, random_state=1),
                 ad_model=Box(), fc_or_rtc_model=ReactionTypeControl()):
        self.metric = metric # возможна только для 'rmse' и для регрессионных метрик только. Невозможно для классификационных метрик
        self.reg_model = reg_model # необходима для нахождения отсечки, если пользователь не задает регресиооную модель,
        # то по умолчанию будет RandomForestRegressor(n_estimators=500, random_state=1). Можно также добавить внутри
        # gridSearchCV если надо
        self.ad_model = ad_model # # необходима для нахождения отсечки, внутреняя модель AD
        self.fc_or_rtc_model = fc_or_rtc_model # необходима для нахождения отсечки, если пользователь не задает AD модель внешнюю

    def fit(self, X, y):
        self.X = check_array(X)
        if y is not None:
            self.y = check_array(y, accept_sparse='csc', ensure_2d=False, dtype=None)  # нужен  для обучения, и
            # для нахождения отсечки
        Y_predict_rf_int = []
        AD_bb_int = []
        for train_index, test_index in rkf.split(X):
            x_train = safe_indexing(X, train_index)
            x_test = safe_indexing(X, test_index)
            y_train = safe_indexing(y, train_index)
            self.clf_int = self.reg_model.fit(x_train, y_train)
            Y_pred_internal = self.clf_int.predict(x_test)
            self.BB_int = self.ad_model.fit(x_train)
            AD_bb_internal = self.BB_int.predict(x_test)
            Y_predict_rf_int.append(Series(Y_pred_internal, index=test_index))
            AD_bb_int.append(Series(AD_bb_internal, index=test_index))
        Y_predict_rf_int = concat(Y_predict_rf_int, axis=1)
        Y_pred_mean_int = np.nanmean(np.array(Y_predict_rf_int), axis=1)
        Y_pred_var_int = np.sqrt(np.nanvar(np.array(Y_predict_rf_int), axis=1))
        AD_bb_int = concat(AD_bb_int, axis=1)
        AD_bb_mean_int = np.nanmean(np.array(AD_bb_int), axis=1)
        self.param_consens = UniqueZ(Y_ext=y, Y_pred=Y_pred_mean_int, AD=AD_bb_mean_int, metric=self.metric, more='>=')
        self.param_var = UniqueZ(Y_ext=y, Y_pred=Y_pred_mean_int, AD=Y_pred_var_int, metric=self.metric, more='<=')
        return self

    def predict(self, X):
        self.X = check_array(X)
        check_is_fitted(self, ['param_consens', 'param_var'])
        Y_predict_rf_ext = []
        AD_bb_ext = []
        for train_index, test_index in rkf.split(X): # надо проверить она будет 50 раз предсказывать свойство?
            Y_pred_external = self.clf_int.predict(X)
            AD_bb_external = self.BB_int.predict(X)
            AD_BB_FC = np.logical_and(AD_bb_external, self.fc_or_rtc_model)
            Y_predict_rf_ext.append(Y_pred_external)
            AD_bb_ext.append(AD_BB_FC)
        Y_pred_var_ext = np.sqrt(np.nanvar(np.array(Y_predict_rf_ext), axis=0))
        AD_bb_mean_ext = np.mean(np.array(AD_bb_ext), axis=0)
        AD_consensus_ext = AD_bb_mean_ext >= self.param_consens['z']
        AD_pred_var = Y_pred_var_ext <= self.param_var['z']
        return AD_consensus_ext, AD_pred_var
