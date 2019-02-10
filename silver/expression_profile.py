import pandas as pd
import numpy as np
import copy


class ExpressionProfile(object):
    """A container class for representing an expression profile.

    """
    def __init__(self, profile):
        assert isinstance(profile, (pd.DataFrame, pd.core.frame.DataFrame))
        self.__profile = profile

    @property
    def shape(self):
        return self.profile.shape

    def __len__(self):
        return self.profile.shape[0]

    def keys(self):
        return iter(self.profile.index.values)

    def items(self):
        for idx, row in self.profile.iterrows():
            yield idx, list(row)

    def values(self):
        for _, row in self.profile.iterrows():
            yield list(row)

    def __getitem__(self, i):
        if not isinstance(i, int):
            raise TypeError('i must be of type integer.')
        if i < 0 or i >= len(self):
            raise IndexError('Index out of bound error.')
        return list(self.profile.iloc[i, :])

    @property
    def columns(self):
        return list(self.profile.columns)

    @property
    def profile(self):
        return self.__profile

    def get(self, identifier):
        if identifier not in self:
            raise KeyError('identifier cannot be found in he Profile object.')
        return list(self.profile.loc[identifier, :])

    def set(self, identifier, expression):
        if identifier not in self:
            raise KeyError('identifier cannot be found in he Profile object.')
        expres = self.get(identifier)
        assert len(expres) == len(expression)
        # for x, y in zip(expres, expression):
        #     if type(x) != type(y):
        #         raise TypeError('expression components must have the same ' +
        #                         'type as expression values in the profile' +
        #                         'for {}.'.format(identifier))
        self.profile.loc[identifier, :] = copy.copy(expression)

    def __contains__(self, identifier):
        return identifier in self.profile.index.values

    def samples(self, cols, is_index=True):
        if is_index is False:
            if not set(cols).issubset(set(self.columns)):
                print(cols, self.columns)
                raise ValueError('cols must be a list of column names.')
            result = self.profile.loc[:, cols]
        else:
            if max(cols) >= len(self.columns) or min(cols) < 0:
                raise IndexError('Index out of bound error.')
            result = self.profile.iloc[:, cols]
        return ExpressionProfile(copy.copy(result))

    def data(self):
        return copy.copy(self.profile.values)

    def __str__(self):
        return str(self.profile)

    def concat(self, other):
        msg = ('Profiles with unequal (different or unaligned) row names'
               'cannot be concatenated')
        assert list(self.keys()) == list(other.keys()), msg
        data = pd.concat([self.profile, other.profile], axis=1)
        return ExpressionProfile(data)