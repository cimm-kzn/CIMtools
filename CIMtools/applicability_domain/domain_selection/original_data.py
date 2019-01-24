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
from CIMtools.conditions_container import ConditionsToDataFrame
from CIMtools.preprocessing import EquationTransformer, SolventVectorizer, Fragmentor, CGR
from sklearn.compose import ColumnTransformer
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from CIMtools.conditions_container import DictToConditions


def fordate(data, constant_name):
    """
    получаем метаданные, откуда берем Y и условия реакции (conditions)
    DictToConditions - это помощник для маппирования словаря условий в объекты условий. Он так же самотестируется,
    то есть он проверит все условия на валидность.
    В DictToConditions реализованы: температура, давление, растворитель и доля растворителя. Для того чтобы им
    воспользоваться, нужно прописать как каждый тип указан в исходных данных, причем в случае растворителей нужно
    указывать лист названий растворителя и названий его доли. Причем, названия растворителей и долей нужно указывать в
    одном порядке. Например, solvents=['additive.1', 'additive.2'] и amounts=['amount.1', 'amount.2'].
    На выходе получаем список объектов условий.

    :param data: data = RDFread('da.rdf').read()
    :param constant_name:
    :return: Y, список структур и условий. В итоге получается матрица длиной N объектов с двумя колонками: структура и
    условие.
    """
    meta = [x.meta for x in data]
    Y = [float(x[constant_name]) for x in meta]
    cond = DictToConditions('temperature', solvents=['additive.1']).transform(meta)
    struct_cond = list(zip(data, cond))
    return Y, struct_cond


def initial_data(data_train, constant_name_train, data_test, constant_name_test): # функция для получения X, y
    """
    ----------------------------- Первый вариант -------------------------------
    Подходит и для того чтобы найти результаты AD внутри одного dataset. Так как для оптимизации гиперпараметров АД
    используем двукратно вложенную валидацию необходимо N раз генерировать X (где N - количество фолдов).
    Пример как выглядит структура кода для одного dataset
    ------------------------------------------------------
    импортируем все необходимые библиотеки и модули. В моем случае так

    path = '/home/asima/asima/FR/datasets' # путь до датасетов
    name = 'SN2_article_2018_final++.rdf' # название датасета
    constant_name = 'logK'
    file = '{}/{}'.format(path, name)
    environ["PATH"]+=":/home/asima" # путь до фрагментора
    for score in ['ba', 'rmse']: # метрика для оптимизации гиперпараметров АД
        AD_gpr = [[] for _ in range(3)] # для gpr нужен свой список так как предсказанные значения Y у GPR свои
        if score == 'ba':
            AD_res = [[] for _ in range(36)] # для каждого AD будет сохраняться значения AD в виде True/False.
            # чтобы потом посчитать общую таблицу
        else:
            AD_res = [[] for _ in range(40)] # 4  метода считаются только при rmse (в самом конце)
        Y = []  # в этом листе сохраним все значения у для расчета полной таблицы
        Y_gpr = []  # тоже самое только для гаусовых процессов так как это отдельная регрисионная модель
        seed = 1  # for random state
        kf = KFold(n_splits=5, random_state=seed, shuffle=True)  # for GridSearchCV
        rkf = RepeatedKFold(n_splits=5, random_state=seed, n_repeats=10)  # for rkf_method (it is one type of AD)
        reactions = RDFread(file).read()

        for train_index_ext, test_index_ext in kf.split(reactions):  # здесь делим датасет на внешнюю обучающую и
        # внешнюю тестовую выборки
            reactions_train = safe_indexing(reactions, train_index_ext) # вненшняя обучающая выборка
            reactions_test = safe_indexing(reactions, test_index_ext) # внешняя тестовая выборка
            data = initial_data(reactions_train, constant_name, reactions_test, constant_name)
            X_train = data[0][0] # np.array
            Y_train = data[0][1] # list
            y_train_gpr = data[0][2][:, 0] # np.array
            X_test = data[1][0]  # np.array
            Y_test = data[1][1]  # list
            y_test_gpr = data[1][2][:, 0]  # np.array

            # Далее строим регресионные модели и АД модели. Надо помнить что для каждого фолда будет свой размер X_train
            # Внутри одного фолда размер X_train.shape == X_test.shape (out: True)
            # Но между фолдами этого может не наблюдаться, так как каждый раз в reactions_train и reactions_test попадают
            # разное количество реакций

    ---------------------------- Второй вариант ---------------------------------------
    Данная функция подходит и для сравнения двух разных датасетов
    Главное на обучающем датасете зафитить фрагментор
    Тестовую выборку отрансформить
    Это нужно для того чтобы проверить сколько реакций чужого типа пропускает AD models
    Пример
    -------
    # import библиотек

    path = '/home/asima/asima/FR/datasets' # путь до датасетов
    name_1 = 'E2_28_08_2018.rdf' # название датасета 1 - обучающая выборка  Нужно будет еще для названия пиклов
    name_2 = 'Tautomres_new_739.rdf' # 2 - тестовая выборка. Нужно будет еще для названия пиклов
    constant_name_train = 'logK' # Таутомеры - 'tabulated_constant'
    constant_name_test = 'tabulated_constant'

    file_train = '{}/{}'.format(path, name_1)
    file_test = '{}/{}'.format(path, name_2)

    environ["PATH"]+=":/home/asima" # путь до фрагментора

    seed = 1 # for random state
    kf = KFold(n_splits=5, random_state=seed, shuffle=True) # for GridSearchCV
    rkf = RepeatedKFold(n_splits=5, random_state=seed, n_repeats=10) # for rkf_method (it is one type of AD)

    reactions_train = RDFread(file_train).read()
    reactions_test = RDFread(file_test).read()
    data = initial_data(reactions_train, constant_name_train, reactions_test, constant_name_test)
    X_train = data[0][0]
    Y_train = data[0][1]
    y_train_gpr = data[0][2][:, 0] # np.array

    X_test = data[1][0]  # np.array
    Y_test = data[1][1]  # list
    y_test_gpr = data[1][2][:, 0]  # np.array

    :param data_train: после того как прочитала реакции reactions = RDFread(file).read(), подаем reactions_train -
                       в первом случае
                       во втором случае надо вначале прочитать reactions_train = RDFread(file_train).read() потом подать
    :param constant_name_train: предсказываемое свойство
    :param data_test: после того как прочитала реакции reactions = RDFread(file).read(), подаем reactions_test -
                       в первом случае
                       во втором случае надо вначале прочитать reactions_test = RDFread(file_test).read() потом подать
    :param constant_name_test: предсказываемое свойство (может совпадать с constant_name_train в первом случае,
                               однако во втором случае могут различаться)
    :return: два list. В первом X, y, y_gpr для обучающей выборки, во втором X, y, y_gpr для тестовой выборки
    """

    data_train = fordate(data_train, constant_name_train)  # из Рамиля workflow для генерации данных
    data_test = fordate(data_test, constant_name_test)
    """
    Далее мы создаем два Pipeline для обработки структур и условий. ColumnTransformer задает определенные трансформеры 
    к каждой колонке. Так как у нас всего две колонки - структуры и условия, то обработка структур будет выглядеть 
    следующим образом: ('structure', Pipeline(...), [0]) где [0] - первая колонка в struct_cond. Аналогично прописываем 
    для условий, только уже указываем индекс для второй колонки. ConditionsToDataFrame - специальный transformer из 
    CIMtools, который превращает список условий в DataFrame, то есть он распаковывает условия в колонки. 
    По умолчанию он вытащит температуру, давление и растворитель. Затем, данный DataFrame преобразуется через 
    следующие трансформеры: EquationTransformer (преобразует температуру к необходимому виду) и SolventVectorizer 
    (превращает название растворителя в 13 дескрипторов). Если мы у нас есть дополнительные дескрипторы доли и 
    растворителя, то можно не использовать transformers, а указать 'passthrough'. Например, для долей растворителя это 
    будет выглядеть так: ('amount', 'passthrough', ['solvent_amount.1'])
    """
    dsc = ColumnTransformer([('conditions', Pipeline([('cd', ConditionsToDataFrame()),
                                                      ('ct', ColumnTransformer([('temperature',
                                                                                 EquationTransformer('1/x'),
                                                                                 'temperature'),
                                                                                ('solvent', SolventVectorizer(),
                                                                                 'solvent.1')]))]), [1]),
                             ('structure', Pipeline([('cgr', CGR()),
                                                     ('frg', Fragmentor(version='2017.x', max_length=4,
                                                                        useformalcharge=True)),
                                                     ('scaler', StandardScaler())]), [0])])  # вначале растворители,
    # потом дескрипторы, дескрипторы в свою очередь стандартизируются на нулевое среднее и единичную дисперсию
    X_train = dsc.fit_transform(data_train[1])
    X_test = dsc.transform(data_test[1])
    # Y for GRP needs scaler
    scaler_gpr = StandardScaler()
    y_train_gpr = scaler_gpr.fit_transform(np.array(data_train[0]).reshape(-1, 1))
    y_test_gpr = scaler_gpr.transform(np.array(data_test[0]).reshape(-1, 1))
    return [X_train, data_train[0], y_train_gpr], [X_test, data_test[0], y_test_gpr]
