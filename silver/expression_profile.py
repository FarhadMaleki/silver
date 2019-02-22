"""This module aims at representing a expression profiles.
"""
import copy
import pandas as pd
import numbers


class ExpressionProfile(object):
    """A class for wrapping gene expression profile data.

    """
    def __init__(self, profile):
        """Initialize an ExpressionProfile object using a pandas DataFrame.

        Args:
            profile(pandas.DataFrame): A DataFrame of expression values.
                Each row represents expression levels for a gene, and each
                column represents expression levels for a sample.

        """
        assert isinstance(profile, (pd.DataFrame, pd.core.frame.DataFrame))
        self.__profile = profile

    @property
    def shape(self):
        """Return the shape of the ExpressionProfile object as a tuple.

        Returns:
            tuple: A tuple of size two, where the first element represents the
                number of genes and the second element represents the number
                of samples.

        """
        return self.profile.shape

    def __len__(self):
        """Return the number of genes in an ExpressionProfile.

        Returns:
            int: The number of genes in the ExpressionProfile.

        """
        return self.profile.shape[0]

    def keys(self):
        """Return the gene IDs (names).

        Returns:
            An iterator of the gene ID/names in the ExpressionProfile object.

        """
        return iter(self.profile.index.values)

    def items(self):
        """Yield tuples of (gene name/ID, expression values).

        Yields:
            tuple: A tuple of size two, where the first element is a
                gene name/ID and the second element is the list of
                expression values, of all samples, for that gene.

        """
        for idx, row in self.profile.iterrows():
            yield idx, list(row)

    def values(self):
        """Yield list of expression values per gene.

        Yields:
            list: the list of expression values, of all samples, for the
                next gene.

        """
        for _, row in self.profile.iterrows():
            yield list(row)

    def __getitem__(self, i):
        """Implement evaluation of self[key].

        Args:
            int: An non-negative number between 0 and self.len()-1
                (both inclusive).

        Returns:
            list: Expression values, of all samples, for the i-th gene.

        """
        if not isinstance(i, numbers.Integral):
            raise TypeError('i must be an integer number.')
        if i < 0 or i >= len(self):
            raise IndexError('Index out of bound error.')
        return list(self.profile.iloc[i, :])

    @property
    def columns(self):
        """ Return the list of sample names.

        Returns:
            list: List of all samples in the ExpressionProfile object.

        """
        return list(self.profile.columns)

    @property
    def profile(self):
        """Return expression profile object as a pandas.DataFrame.

        Args:
            pandas.DataFrame: Expression values of all genes across all samples.
        """
        return self.__profile

    def get(self, identifier):
        """Get the expression values for a gene.

        Args:
            str: A gene name/ID.

        Returns:
            list: A list of expression values for the gene with identifier as
                its name/ID

        Raises:
            KeyError: Raise KeyError if identifier is not the name/ID of a gene
                in the ExpressionProfile object.

        """
        if identifier not in self:
            raise KeyError('identifier cannot be found in the Profile object.')
        return list(self.profile.loc[identifier, :])

    def set(self, identifier, expression):
        """Set the expression values of a gene.

        Args:
            identifier (str): A gene name/ID
            expression (list): A list of expression values to be assigned to
                the gene with identifier as its name. The length of the
                expression argument must be equal to the number of samples in
                the ExpressionProfile.

        Raises:
            KeyError: Raise KeyError if identifier is not the name/ID of a gene
                in the ExpressionProfile object.
            ValueError: Raise ValueError if the length of expression is not
                equal to the number of samples in the ExpressionProfile.

        """
        if identifier not in self:
            raise KeyError('identifier cannot be found in the Profile object.')
        expres = self.get(identifier)
        if len(expres) != len(expression):
            raise ValueError('length of expression must be equal to the ' +
                             'number of samples in the ExpressionProfile.')
        self.profile.loc[identifier, :] = copy.copy(expression)

    def __contains__(self, identifier):
        """Check if a gene exist in the ExpressionProfile object.

        Args:
            identifier (str): A gene name/ID.

        Returns:
            bool: True if identifier is a gene name/ID in the ExpressionProfile
                object; False otherwise.

        """
        return identifier in self.profile.index.values

    def samples(self, cols):
        """Create an ExpressionProfile from a subset of samples.

        Args:
            cols (list): A sequence of column names or column indices.
                If cols is a list of column names/IDs, it must be a subset of
                sample names/IDs in the ExpressionProfile. If cols is a list of
                integers its elements must be between 0 and  the number of
                samples -1.

        Returns:
            ExpressionProfile: An ExpressionProfile object containing the
                expression level of all genes for the samples represented
                in cols.

        Raises:
            ValueError: Raises ValueError if cols is a list of the
                ExpressionProfile sample  names/IDs (of type strings) or
                a list of integers between 0 and  the number of samples -1.
            IndexError: Raises IndexError if cols contains indices that are
                out of bound, i.e. not between 0 and  the number of samples -1.

        """
        error_msg = ('cols must be a list of the ExpressionProfile sample '
                     'names/IDs (of type strings) or a list of integers '
                     'between 0 and  the number of samples -1.')
        is_index = isinstance(cols[0], numbers.Integral)
        #All elements of cols must have the same data type
        if all(isinstance(col, type(cols[0])) for col in cols) is not True:
            raise ValueError(error_msg)
        if is_index is False:
            if not set(cols).issubset(set(self.columns)):
                raise ValueError(error_msg)
            result = self.profile.loc[:, cols]
        else:
            if max(cols) >= len(self.columns) or min(cols) < 0:
                raise IndexError('Index out of bound error.')
            result = self.profile.iloc[:, cols]
        return ExpressionProfile(copy.copy(result))

    def data(self):
        """Return the expression values of the ExpressionProfile.

        Returns:
            numpy.ndarray: A shallow copy of the expression values from
                the ExpressionProfile object.

        """
        return copy.copy(self.profile.values)

    def __str__(self):
        """A string representation.

        Returns:
            str: A string representation of the ExpressionProfile object.
        """
        return str(self.profile)

    def concat(self, other):
        """Concatenate two ExpressionProfiles.

        Args:
            other (ExpressionProfile): An ExpressionProfile with the same
                keys() and with different sample names comparing to self.

        Returns:
            ExpressionProfile: An ExpressionProfile resulting from
            concatenation of self and other.

        Raises:
            ValueError: Raise ValueError if self and other have different
                or unaligned keys(), i.e. gene names/IDS or if the sample
                names in self and other are not different from each other.

        """
        msg = ('Profiles with unequal (different or unaligned) row names '
               'cannot be concatenated.')
        if  list(self.keys()) != list(other.keys()):
            raise ValueError(msg)
        if set(self.columns) & set(other.columns) != set():
            raise ValueError('Profile sample names must be unique after '
                             'concatenation.')
        data = pd.concat([self.profile, other.profile], axis=1)
        return ExpressionProfile(data)
