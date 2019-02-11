"""Represent a repository, which is a collection for expression values.
"""
import numpy as np


class Repository(object):
    """ Create a repository using case sample expression profile.

    The generated repository contains num_sim_cases and
    num_repetitions * len(case_expression_profile) rows.

        Args:
            case_expression_profile (ExpressionProfile): Expression data for
                case samples.
            num_sim_cases (int): A positive integer representing the number of
                control samples to be simulated.
            num_repetitions (int): A positive integer representing the number of
                samplings. In each sampling a set of columns from the cases is
                extracted and appended to the repository.

        Raises:
            ValueError: If num_sim_cases is greater than the number of samples
                in case_expression_profile.

    """

    def __init__(self, case_expression_profile, num_sim_cases,
                 num_repetitions=1, random_state=None):
        np.random.seed(random_state)
        total_num_cases = len(case_expression_profile.columns)
        if num_sim_cases > total_num_cases:
            raise ValueError('num_sim_cases must be less than or equal to' +
                             ' the total number of case samples.')
        indices = np.random.choice(np.arange(total_num_cases), num_sim_cases,
                                   replace=False)
        indices = [int(value) for value in indices]
        repository = case_expression_profile.samples(indices).data()
        for _ in range(1, num_repetitions):
            indices = list(np.random.choice(np.arange(total_num_cases),
                                            num_sim_cases, replace=False))
            indices = [int(value) for value in indices]
            batch = case_expression_profile.samples(indices).data()
            repository = np.concatenate([repository, batch], axis=0)
        self.repository = repository

    @property
    def shape(self):
        """Return the shape of the Repository object.

        Returns:
            tuple: A tuple of length two, where the first element is
                the number of rows in the Repository object and the second
                element is the number of samples.
        """
        return self.repository.shape

    def __len__(self):
        """Return the length of the repository object.

        Args:
            int: The number of rows in the Repository object.

        """
        return self.repository.shape[0]

    def __getitem__(self, i):
        """Return a sequence of expression values.

        Args:
            i (int): The index of a row of the Repository object.
                It must be between 0 and len(self)-1 (both inclusive).

        Returns:
            array-like:Return the sequence of expression values in the
                i-th row of the Repository object.

        Raises:
            IndexError: If i is not between 0 and len(self)-1
                (both inclusive).
        """
        if 0 <= i < len(self):
            return self.repository[i, :]
        else:
            raise IndexError('Index out of bound.')

    def __iter__(self):
        """Create an iterator of the repository rows.

        Returns:
            iterator: An iterator of the Repository rows.

        """
        return iter(self.repository)

    def __str__(self):
        """A String reporesentation.

        Return:
            str: A string represetntaion of the Repository object.

        """
        return str(self.repository)
