from CGRtools.preparer import CGRpreparer
from collections import defaultdict
from numpy import array
from sklearn.model_selection import BaseCrossValidator
from sklearn.utils.validation import indexable


class LeaveOneGroupOut(BaseCrossValidator):
    """Leave-One-Group-Out cross-validator

    Provides train/test indices to split data in train/test sets. Each
    reactions with the same condition is used once as a test set (singleton)
    while the remaining reactions form the training set. Test set includes
    only reactions with transformations that appeared in other reactions.
    """

    def __init__(self):
        # We need this for the build_repr to work properly in py2.7
        # see #6304
        pass

    def get_n_splits(self, X=None, y=None, groups=None):
        """Returns the number of splitting iterations in the cross-validator
        Parameters
        ----------
        X : object
            Always ignored, exists for compatibility.
            ``np.zeros(n_samples)`` may be used as a placeholder.
        y : object
            Always ignored, exists for compatibility.
            ``np.zeros(n_samples)`` may be used as a placeholder.
        groups : array-like, with shape (n_samples,), optional
            Group labels for the samples used while splitting the dataset into
            train/test set.
        Returns
        -------
        n_splits : int
            Returns the number of splitting iterations in the cross-validator.
        """
        return len(set(groups))

    def split(self, X, y=None, groups=None):
        """Generate indices to split data into training and test set.
        Parameters
        ----------
        X : array-like, of length n_samples
            Training data, includes reaction's containers
        y : array-like, of length n_samples
            The target variable for supervised learning problems.
        groups : array-like, with shape (n_samples,), optional
            Group labels for the samples used while splitting the dataset into
            train/test set.
        Yields
        ------
        train : ndarray
            The training set indices for that split.
        test : ndarray
            The testing set indices for that split.
        """
        X, y, groups = indexable(X, y, groups)
        cgr = CGRpreparer()
        cgrs = [cgr.condense(r) for r in X]

        structure_condition = defaultdict(set)

        for structure, condition in zip(cgrs, groups):
            structure_condition[structure].add(condition)

        train_data = defaultdict(list)
        test_data = []

        for n, (structure, condition) in enumerate(zip(cgrs, groups)):
            train_data[condition].append(n)
            if len(structure_condition[structure]) > 1:
                test_data.append(n)

        for condition, indexes in train_data.items():

            test_index = [index for index in indexes if index in test_data]
            if len(test_index) == 0:
                continue

            train_index = []
            for c in train_data.keys():
                if not c == condition:
                    train_index.extend(train_data[c])

            yield array(train_index), array(test_index)
