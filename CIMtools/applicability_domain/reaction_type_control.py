# -*- coding: utf-8 -*-
#
#
#  Copyright 2018 Assima Rakhimbekova <asima.astana@outlook.com>
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


class ReactionTypeControl():
    """
    создаем класс, который создает список хэш-ключей реакций с задаваемым нами радиусом сферы
    этот список и является областью применимости модели, то есть фитинг - это просто список ключей, а проверяет есть ли
    такие фрагменты в домейне или нет, на выходе дает ответ True or False.

    NB! по умолчанию у меня гибридизация True! Если адо будет перебирать его как гиперпараметр надо изменить код
    """
    def __init__(self, env=0, hyb=False):
        self.env = env
        self.hyb = hyb

    def __readFile(self, X, index=None):
        if index is not None:  # бывают случаю когда надо фитить по фолдам, тогда необходимо это
            if self.env == 'all':
                return [str(~r) for num, r in enumerate(X) if num in index]
            else:
                data = []
                for num, r in enumerate(X):
                    if num in index:
                        cgr = ~r
                        cgr.reset_query_marks()
                        center_atoms = cgr.get_center_atoms()
                        aug_center = cgr.augmented_substructure(center_atoms, dante=False, deep=self.env,
                                                                as_view=True)
                        data.append("{:h}".format(aug_center))
        else:
            if self.env == 'all':
                return [str(~r) for r in X]
            else:
                data = []
                for r in X:
                    cgr = ~r  # Condence Graph of Reaction
                    cgr.reset_query_marks()  # reset hyb and neighbors marks to atoms
                    center_atoms = cgr.get_center_atoms()  # Numbers atoms in reaction center in list
                    aug_center = cgr.augmented_substructure(center_atoms, dante=False, deep=self.env,
                                                            as_view=True)  # get ubgraph with atoms and their neighbors
                    data.append("{:h}".format(aug_center))  # String for graph reaction center

        return data

    def fit(self, X, index=None):
        """Fit distance-based AD.

        Parameters
        ----------
        X : rdf file.
        index : в случае если надо будет фитить по фолдам
        Returns
        -------
        self : object
            Returns self.
        """
        dataHash = set(self.__readFile(X, index))
        self.bytes = [bytes(list(dataHash)[i], encoding='utf-8') for i in range(len(dataHash))]
        return self

    def predict(self, X, index=None):
        """

        :param X: rdf file.
        :param index: в случае если надо будет предсказывать по фолдам
        :return: результат объекты принадлежат области применимости модели или нет. Если хеш объекта из тестовой выборки
        совпадает со списком уникальных хэшей из обучающей выборки, то объект принадлежит области применимости модели, в
        противном случае нет.
        """
        state = [bytes(newHash, encoding='utf-8') in self.bytes for newHash in self.__readFile(X, index)]
        return state


__all__ = [ReactionTypeControl]
