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
import pandas as pd
from sklearn.metrics import mean_squared_error, confusion_matrix, roc_curve, auc, r2_score
from sklearn.model_selection import KFold, RepeatedKFold


seed = 1
kf = KFold(n_splits=5, random_state=seed, shuffle=True)
rkf = RepeatedKFold(n_splits=5, random_state=seed, n_repeats=10)

def res(AD_list, num, AD, FC=None, RTC_deep_1=None):
    """
    Эта функция добавляет во внешний лист значения AD чтобы потом посчитать суммарную таблицу
    path = '/home/asima/asima/FR/datasets' # путь до датасетов
    name = 'SN2_article_2018_final++.rdf' # название датасета 1 - обучающая выборка  Нужно будет еще для названия пиклов
    constant_name = 'logK' # Таутомеры - 'tabulated_constant'
    file = '{}/{}'.format(path, name)
    environ["PATH"]+=":/home/asima" # путь до фрагментора
    for score in ['ba', 'rmse']:
        print(score)
        AD_gpr = [[] for _ in range(3)] # для gpr нужен свой список так как предсказанные значения Y у GPR свои
        if score == 'ba':
            AD_res = [[] for _ in range(36)] # для каждого AD будет сохраняться значения AD в виде True/False. чтобы потом посчитать общую таблицу
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
        ... строим модель ...
        # **************************************** REACTION TYPE CONTROL DEEP = 1 **************************************
        RTC_1 = ReactionTypeControl(env=1, hyb=True).fit(reactions_train)
        AD_RTC_deep1 = RTC_1.predict(reactions_test)
        res(AD_res, 0, AD_RTC_deep1) # В list AD_res add results AD_RTC_deep1 in AD_res[0]
        -------
        # *********************************** RFR_VAR ******************************************************************
        est_var = RandomForestRegressor2(random_state=seed, n_estimators=500, max_features=est.best_params_['max_features'],
                                         n_jobs=4).fit(X_train, Y_train)
        AD_est_var_values = est_var.predict_proba(X_test)
        min_h_param_RFR_VAR = threshold(ad='rfr_var', X=X_train, y=Y_pr_ts, metric=score, reg_model=est.best_estimator_,
                                        ad_model=est_var) # для нахождения отсечки
        AD_var_RF_itog = AD_est_var_values <= min_h_param_RFR_VAR['z']
        res(AD_res, 2, AD_var_RF_itog, FC=AD_FC, RTC_deep_1=AD_RTC_deep1)
    :param AD_list: название листа куда будут сохраняться данные
    :param num: в какой лист будет записываться
    :param AD: результаты АД
    :param FC: результат FC
    :param RTC_deep_1:  results RTC_deep_1
    :return: заполненный лист со всеми адшками. Причем обязательно надо extend-ить так как если делим на две кросс-валидации
    надо сохранить в том же порядке что и у предсказзанные, у_true and ad. Все должно быть идентично
    """
    AD_list[num].extend(AD)
    if FC is not None:
        AD_with_FC = np.logical_and(AD, FC)
        AD_list[num + 1].extend(AD_with_FC) # комбинированный метод с FC Если и AD и RTC_deee_1 говорят что объект в AD, то объект в AD
    if RTC_deep_1 is not None:
        AD_with_RTC = np.logical_and(AD, RTC_deep_1) # комбинированный метод с RTC. Если и AD и RTC_deee_1 говорят что
        # объект в AD, то объект в AD
        AD_list[num + 2].extend(AD_with_RTC)
    print('is_done!')


def my_metrics_(Y_exp, Y_pred, AD_pred):
    """
    Метрики для характеризации качества предсказательной модели
    для определения лучшего АД
    :param Y_exp: экспериментальные значения у (из контрольной внешней выборки)
    :param Y_pred: предсказанные значения у случайным лесом
    :param AD_pred: АД
    :return: R2_noAD, R2_AD, R2_AD_out_n, R2_score, R2_score_out_n, diff_R2, RMSE_noAD, RMSE_AD, RMSE_AD_out_n, \
           RMSE_score, RMSE_score_out_n, diff_RMSE, coverage
           R2_noAD - R2 без учета АД
           R2_AD - с учетом АД то есть с учетом только инлаеров
           R2_AD_out_n - с учетом outlier
           R2_score - разница R2_noAD-R2_AD
           и другие метрики ...
    """
    AD_pred_n = (AD_pred == 1)
    AD_pred_out = AD_pred_n - 1
    AD_pred_out_n = (AD_pred_out == -1)

    R2_noAD = r2_score(Y_exp, Y_pred)
    average_Y_exp = np.mean(Y_exp)
    numerator = sum(map(lambda x: (((x[0] - x[1]) ** 2) * x[2]), zip(Y_pred, Y_exp, AD_pred_n)))
    denominator = sum([((x[0] - average_Y_exp) ** 2) * x[1] for x in zip(Y_exp, AD_pred_n)])
    numerator_out_n = sum(map(lambda x: (((x[0] - x[1]) ** 2) * x[2]), zip(Y_pred, Y_exp, AD_pred_out_n)))
    denominator_out_n = sum([((x[0] - average_Y_exp) ** 2) * x[1] for x in zip(Y_exp, AD_pred_out_n)])
    if numerator and denominator:
        R2_AD = 1 - (numerator / denominator)
    else:
        R2_AD = 0
    R2_score = R2_AD - R2_noAD

    if numerator_out_n and denominator_out_n:
        R2_AD_out_n = 1 - (numerator_out_n / denominator_out_n)
    else:
        R2_AD_out_n = 0
    R2_score_out_n = R2_AD_out_n - R2_noAD
    diff_R2 = R2_AD_out_n - R2_AD

    RMSE_noAD = np.sqrt(mean_squared_error(Y_exp, Y_pred))
    # RMSE_noAD = sqrt((sum(map(lambda x: (x[0] - x[1]) ** 2, zip(Y_pred, Y_exp)))) / len(Y_exp))
    s_n = sum(AD_pred_n)
    s_out_n = sum(AD_pred_out_n)
    if s_n:
        RMSE_AD = np.sqrt((sum(map(lambda x: (((x[0] - x[1]) ** 2) * x[2]), zip(Y_pred, Y_exp, AD_pred_n)))) / s_n)
    else:
        RMSE_AD = 0
    RMSE_score = RMSE_AD - RMSE_noAD

    if s_out_n:
        RMSE_AD_out_n = np.sqrt(
            (sum(map(lambda x: (((x[0] - x[1]) ** 2) * x[2]), zip(Y_pred, Y_exp, AD_pred_out_n)))) / s_out_n)
    else:
        RMSE_AD_out_n = 0
    RMSE_score_out_n = RMSE_AD_out_n - RMSE_noAD
    diff_RMSE = RMSE_AD - RMSE_AD_out_n

    coverage = sum(AD_pred_n) / len(Y_exp)
    return R2_noAD, R2_AD, R2_AD_out_n, R2_score, R2_score_out_n, diff_R2, RMSE_noAD, RMSE_AD, RMSE_AD_out_n, \
           RMSE_score, RMSE_score_out_n, diff_RMSE, coverage


def my_confusion_matrix(Y_pred, AD_pred):
    """
    Метрика №2
    :param Y_pred: Y_pred = abs(Y_pred_test[:, 0] - Y_pred_test[:, 1]) <= 3 * np.sqrt(mean_squared_error(Y_pred_test[:, 1],
                                                                                          Y_pred_test[:, 0])) то есть значения true/false
                                                                                          тип как "истинные значения"
    :param AD_pred: предсказанные значения true/false (то есть результаты ад)
    :return: tn, fp, fn, tp, acc, sen_TPR, spc_TNR, PPV, PNV, F_score_N, F_score_P, ba
    """
    AD_pred_n = (AD_pred == 1)
    tn, fp, fn, tp = confusion_matrix(Y_pred, AD_pred_n).ravel()
    acc = (float(tp + tn)) / (tn + fp + fn + tp)
    sen_TPR = (float(tp)) / (tp + fn)
    spc_TNR = (float(tn) / (tn + fp))
    ba = ((sen_TPR + spc_TNR) / 2)
    PPV = (float(tp)) / (tp + fp)
    PNV = (float(tn)) / (tn + fn)
    F_P = 1 / PPV + 1 / sen_TPR
    F_score_P = 1 / F_P
    F_N = 1 / PNV + 1 / spc_TNR
    F_score_N = 1 / F_N
    return tn, fp, fn, tp, acc, sen_TPR, spc_TNR, PPV, PNV, F_score_N, F_score_P, ba


def IAP(Y_pred, Y_exp, AD_pred):
    """
    IAP - Invariant Accuracy of Prediction
    Еще одна метрика для описания качества построенной модели
    :param Y_pred: предсказанные значения с помощью случайного леса на тестовом сэте
    :param y_test: истинные значения y_test
    необходимы они для того чтобы посчитать ошибку предсказания значений
    :param AD_pred: True или False то бишь входит в АД реакция из тест сета или нет
    Затем создаем матрицу из N строк и 2х столбцов
    Первый столбец - ошибка предсказания
    Второй столбец - АД
    затем сортируем эту матрицу по возрастанию ошибки, то есть вначале маленькая ошибка а в конце самая большая оишбка
    АД должны меняться синхронно сортировке
    Получается что True больше сверху, а False снизу
    Затем считаем сумму сколько False после кажддой True и делим на произведение сколько True и cколько False
    :return: IAP

    Например:
    ---------
    0 True
    0,1 False
    0,2 True
    0,3 False
    0,4 False
    IAP = (3+2)/(2*3)
    Находим первую True а потом считаем сколько False
    Далее ищем следующую True а потом снова считаем сколько встретяться False и тд
    А затем делим на количество True * количество False
    """
    AD_pred_n = (AD_pred == 1)
    table = np.column_stack((np.array(abs(Y_pred - Y_exp)), AD_pred_n))
    table_sort = table[np.lexsort(np.fliplr(table).T)]  # это какая-то шайтан функция но делает то что я хочу
    only_bool = table_sort[:, 1].tolist()  # оставляем только АД
    if only_bool.count(True) and only_bool.count(False):  # проверяем условие !нам надо чтобы были и True and False!
        i = 0
        b_list = []
        while True:
            try:
                i = only_bool.index(1, i) + 1
                b_list.append(only_bool[i:].count(0))
            except:
                break
        iap = sum(b_list) / (only_bool.count(0) * only_bool.count(1))
    else:
        iap = None
    return iap


def my_auc(AD, Y_probability):
    """
    :param AD: результаты ад в виде true/false
    :param Y_probability: Y_probability = -abs(Y_pred_rf - Y_test)  # FOR AUC
    :return: roc_auc
    """
    AD_pred = (AD == 1)
    fpr, tpr, thresholds = roc_curve(AD_pred, Y_probability)
    roc_auc = auc(fpr, tpr)
    return roc_auc


def result_for_all(i, name_list, Y, AD_pred):
    AD = (AD_pred == 1)
    result_for_fold_1 = my_metrics_(Y_exp=Y[:, 1], Y_pred=Y[:, 0], AD_pred=AD)
    name_list[i].extend(result_for_fold_1)
    try:
        result_for_fold_2 = my_confusion_matrix(Y_pred=Y[:, 2], AD_pred=AD)
        name_list[i].extend(result_for_fold_2)
    except ValueError:
        name_list[i].extend('None')
    result_for_fold_3 = IAP(Y_pred=Y[:, 0], Y_exp=Y[:, 1], AD_pred=AD)
    name_list[i].append(result_for_fold_3)
    result_for_fold_4 = my_auc(AD, Y[:, 3])
    name_list[i].append(result_for_fold_4)
    return name_list


def result(Y, AD, score):
    """
    На выходе будет такая структура результата
	                1R2_noAD	2R2_AD	3R2_AD_out_n 4R2_score 5R2_score_out_n 6diff_R2	7RMSE_noAD 8RMSE_AD	 9RMSE_AD_out_n 10RMSE_score ...
    RTC (deep=1)	0,8057713	0,8152132	0,4701428	0,009442	-0,335628	-0,34507	0,5145257	0,4989964	1,1060841	-0,015529	0,5915584	-0,607088	0,9838476	11	84	67	4667	0,9687306	0,9858471	0,1157895	0,9823195	0,1410256	0,0635838	0,4920401	0,5508183	0,8082617	0,8082617
    FC	0,8057713	0,8124337	0,3679477	0,0066625	-0,437824	-0,444486	0,5145257	0,5042804	1,1518484	-0,010245	0,6373226	-0,647568	0,9902671	7	88	40	4694	0,9734935	0,9915505	0,0736842	0,9815977	0,1489362	0,0492958	0,4932745	0,5326173	0,7476263	0,7476263
    RFR_VAR	0,8057713	0,8063599	0,7634603	0,0005886	-0,042311	-0,0429	0,5145257	0,5107911	1,3935367	-0,003735	0,8790109	-0,882746	0,9977221	2	93	9	4725	0,9788776	0,9980989	0,0210526	0,9806974	0,1818182	0,0188679	0,4946608	0,5095757	0,9071852	0,9071852
    RFR_VAR&FC	0,8057713	0,8131255	0,5569647	0,0073543	-0,248807	-0,256161	0,5145257	0,5004075	1,2014277
    :param Y: все у собранные после всех фолдов (регрисионная модель)
    :param AD: все ад собранные после всех фолдов
    :param score: метрика по которой оптимизировать
    :return: таблица
    """
    Y_np = np.array(Y)
    AD_np = np.array(AD)
    ii = len(AD)
    big = [[] for _ in range(ii)]
    FPR = [[] for _ in range(ii)]
    TPR = [[] for _ in range(ii)]
    for i in range(ii):
        call_all = result_for_all(i=i, name_list=big, Y=Y_np, AD_pred=np.array(AD_np[i]))
        fpr, tpr, _ = roc_curve(np.array(AD[i]), Y_np[:, 3])
        FPR[i].append(fpr)
        TPR[i].append(tpr)
    column_index = np.array(['1R2_noAD', '2R2_AD', '3R2_AD_out_n', '4R2_score', '5R2_score_out_n', '6diff_R2',
                             '7RMSE_noAD', '8RMSE_AD', '9RMSE_AD_out_n', '10RMSE_score', '11RMSE_score_out_n',
                             '12diff_RMSE', '13coverage', '14tn', '15fp', '16fn', '17tp', '18acc', '19sen_TPR',
                             '20spc_TNR', '21PPV', '22PNV', '23F_score_N', '24F_score_P', '25ba', '26IAP', '27auc'])
    if score == 'ba':
        index = ['RTC (deep=1)', 'FC', 'RFR_VAR', 'RFR_VAR&FC', 'RFR_VAR&RTC', 'RTC', 'BB', 'BB&FC', 'BB&RTC', '2CC',
                 '2CC&FC', '2CC&RTC', 'Lev&cv', 'Lev&cv&FC', 'Lev_cv&RTC', 'Leverage', 'Leverage&FC', 'Leverage&RTC',
                 'Z1NN_cv', 'Z1NN_cv&FC', 'Z1NN_cv&RTC', 'Z1NN', 'Z1NN&FC', 'Z1NN&RTC', '1-SVM', '1-SVM&FC',
                 '1-SVM&RTC',
                 'IF', 'IF&FC', 'IF&RTC', 'OZ', 'PZ', 'R', 'IdealM', 'Hybrid_1 (mean)', 'Hybrid_2 (1)']
    else:
        index = ['RTC (deep=1)', 'FC', 'RFR_VAR', 'RFR_VAR&FC', 'RFR_VAR&RTC', 'RTC', 'BB', 'BB&FC', 'BB&RTC', '2CC',
                 '2CC&FC', '2CC&RTC', 'Lev&cv', 'Lev&cv&FC', 'Lev_cv&RTC', 'Leverage', 'Leverage&FC', 'Leverage&RTC',
                 'Z1NN_cv', 'Z1NN_cv&FC', 'Z1NN_cv&RTC', 'Z1NN', 'Z1NN&FC', 'Z1NN&RTC', '1-SVM', '1-SVM&FC',
                 '1-SVM&RTC',
                 'IF', 'IF&FC', 'IF&RTC', 'OZ', 'PZ', 'R', 'IdealM', 'Hybrid_1 (mean)', 'Hybrid_2 (1)', 'Consensus_FC',
                 'PredVar_FC', 'Consensus_RTC', 'PredVar_RTC']
    big_data = pd.DataFrame(big, index=index, columns=column_index)
    big_data.to_excel('{}.xlsx'.format(score))


def result_gpr(Y, AD, score):
    """
    На выходе будет такая структура результата
	                1R2_noAD	2R2_AD	3R2_AD_out_n 4R2_score 5R2_score_out_n 6diff_R2	7RMSE_noAD 8RMSE_AD	 9RMSE_AD_out_n 10RMSE_score ...
    GPR	0,8057713	0,8152132	0,4701428	0,009442	-0,335628	-0,34507	0,5145257	0,4989964	1,1060841	-0,015529	0,5915584	-0,607088	0,9838476	11	84	67	4667	0,9687306	0,9858471	0,1157895	0,9823195	0,1410256	0,0635838	0,4920401	0,5508183	0,8082617	0,8082617
    GPR_FC	0,8057713	0,8124337	0,3679477	0,0066625	-0,437824	-0,444486	0,5145257	0,5042804	1,1518484	-0,010245	0,6373226	-0,647568	0,9902671	7	88	40	4694	0,9734935	0,9915505	0,0736842	0,9815977	0,1489362	0,0492958	0,4932745	0,5326173	0,7476263	0,7476263
    GPR_RTC	0,8057713	0,8063599	0,7634603	0,0005886	-0,042311	-0,0429	0,5145257	0,5107911	1,3935367	-0,003735	0,8790109	-0,882746	0,9977221	2	93	9	4725	0,9788776	0,9980989	0,0210526	0,9806974	0,1818182	0,0188679	0,4946608	0,5095757	0,9071852	0,9071852
    :param Y: все у собранные после всех фолдов
    :param AD: все ад собранные после всех фолдов
    :param score: метрика по которой оптимизировать
    :return: таблица
    """
    Y_np = np.array(Y)
    AD_np = np.array(AD)
    ii = len(AD)
    big = [[] for _ in range(ii)]
    FPR = [[] for _ in range(ii)]
    TPR = [[] for _ in range(ii)]
    for i in range(ii):
        call_all = result_for_all(i=i, name_list=big, Y=Y_np, AD_pred=np.array(AD_np[i]))
        fpr, tpr, _ = roc_curve(np.array(AD[i]), Y_np[:, 3])
        FPR[i].append(fpr)
        TPR[i].append(tpr)
    column_index = np.array(['1R2_noAD', '2R2_AD', '3R2_AD_out_n', '4R2_score', '5R2_score_out_n', '6diff_R2',
                             '7RMSE_noAD', '8RMSE_AD', '9RMSE_AD_out_n', '10RMSE_score', '11RMSE_score_out_n',
                             '12diff_RMSE', '13coverage', '14tn', '15fp', '16fn', '17tp', '18acc', '19sen_TPR',
                             '20spc_TNR', '21PPV', '22PNV', '23F_score_N', '24F_score_P', '25ba', '26IAP', '27auc'])
    index = ['GPR', 'GPR_FC', 'GPR_RTC']
    big_data = pd.DataFrame(big, index=index, columns=column_index)
    big_data.to_excel('GPR{}.xlsx'.format(score))
