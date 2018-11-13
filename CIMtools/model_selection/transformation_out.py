from sklearn.model_selection import BaseCrossValidator
from CGRtools.preparer import CGRpreparer
from collections import defaultdict
from random import shuffle as r_shuffle
from random import random
import numpy as np
from sklearn.utils.validation import indexable


class TransformationOut(BaseCrossValidator):
    def __init__(self, n_splits=5, n_repeats=5, shuffle=False):
        self.n_splits = n_splits
        self.shuffle = shuffle
        self.n_repeats = n_repeats

    def split(self, X, y=None, groups=None):
        X, y, groups = indexable(X, y, groups)
        cgr = CGRpreparer()
        cgrs = [cgr.condense(r) for r in X]

        condition_structure = defaultdict(set)

        for structure, condition in zip(cgrs, groups):
            condition_structure[condition].add(structure)

        train_data = defaultdict(list)
        test_data = []

        for n, (structure, condition) in enumerate(zip(cgrs, groups)):
            train_data[structure].append(n)
            if len(condition_structure[condition]) > 1:
                test_data.append(n)

        if self.n_splits > len(train_data):
            raise ValueError("Cannot have number of splits n_splits=%d greater"
                             " than the number of transformations: %d."
                             % (self.n_splits, len(train_data)))

        structures_weight = sorted({x: len(y) for x, y in train_data.items()}, key=lambda z: z[1], reverse=True)
        fold_mean_size = len(cgrs) // self.n_splits

        folds = [[] for _ in range(self.n_splits)]


        for idx in range(self.n_repeats):
            for i in structures_weight:
                if self.shuffle:
                    r_shuffle(folds)
                for n, fold in enumerate(folds):
                    if len(fold) + i[1] <= fold_mean_size:
                        fold.extend(train_data[i[0]])
                        break
                    else:
                        roulette_param = (len(train_data[i[0]]) - fold_mean_size + len(fold)) / len(train_data[i[0]])
                        if random() > roulette_param:
                            fold.extend(train_data[i[0]])
                            break
                        elif n == 4:
                            fold.extend(train_data[i[0]])


            for i in folds:
                train_index = []
                test_index = []
                for j in range(self.n_splits):
                    if not j == i:
                        train_index.extend(folds[j])
                    else:
                        for c in folds[j]:
                            if groups is not None:
                                if c in test_data:
                                    test_index.append(c)
                            else:
                                test_index.append(c)
                yield np.array(train_index), np.array(test_index)