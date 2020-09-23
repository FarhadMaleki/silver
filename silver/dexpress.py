"""This module contains criteria used for simulating differential expression.

Each criterion must extend DExpress class, which is an abstract class.

"""
from abc import ABC, abstractmethod
import numpy as np
from scipy import stats
from silver.exceptions import NonExistingExpressionException, InvalidLogFoldChangeValue


class DExpress(ABC):
    """An abstract class to be extended by any expression Criterion.

    Args:
        repository:
            Repository: A repository containing the expression values of
                genes for several samples.

    """

    def __init__(self, repository):
        self.repository = repository

    @abstractmethod
    def __call__(self, gene_ctrl_expressions, **kwargs):
        """Simulate differntial expression.

        Args:
            gene_ctrl_expressions(array-like): A sequence of expression values.
            **kwargs: Arguments to be used to simulate differential expression
                of a gene.

        """

class TTestDExpress(DExpress):
    """A criteria based on independent two-sample t-test.

    Args:
        repository (Repository): A expression repository to be used for
            simulating differential expression.
        alpha (float): A positive number between 0 and 1 used as
            significant level.

    """

    def __init__(self, repository, alpha):
        super().__init__(repository)
        self.alpha = alpha

    def __call__(self, gene_ctrl_expressions, fold_change, force=True,
                 shuffle=True, std=0.5, log_scale=True):
        """Simulate differential expression.

        Args:
            gene_ctrl_expressions(array-like): A sequence of expression values.
            fold-change (float): A fold-change value. Note that this value is
                interpreted along with the 'log_scale' parameter. Expression data are
                commonly log normalized to deal with heterogeneity of variance.
                In such cases the fold_change value should be in log scale.
            force (bool): Force the method to return an array of expression
                values. If there is no row in the repository that meet the
                criteria, gene_ctrl_expressions will be shifted using normally
                distributed values where the mean is a random value in the
                fold-change interval and a standard deviation of std.
            shuffle(bool): If True (default), rows of repository will be
                searched in random order for finding the row that meets the
                differential expression criteria; otherwise, the repository
                will be searched from the first row to the last row.
            std(float): A positive value respresenting standard deviation of
                the normally distributed shifts. Only applicable when force
                is True.
            log_scale (bool): True if expression data is in log scale, i.e. data has been
                log normalized, False otherwise.

        Returns:
            numpy.ndarray: An element of the respository that meets the
                differential expression criteria.

        Raises:
            TypeError: If fold_change is not a tuple of length two.
            ValueError: If fold_change[0] >= fold_change[1].

        """
        # Check fold_change type
        if not isinstance(fold_change, tuple) or len(fold_change) != 2:
            raise TypeError('fold_change must be a tuple of length 2.')
        if fold_change[0] > fold_change[1]:
            msg = 'fold_change[0] must be less than fold_change[1]'
            raise InvalidLogFoldChangeValue(msg)
        if fold_change[0] * fold_change[1] <= 0:
            msg = ('fold_change must be a tuple of two numbers of the '
                   'same sign.')
            raise InvalidLogFoldChangeValue(msg)
        if log_scale is True:
            lower = np.array(gene_ctrl_expressions) + fold_change[0]
            upper = np.array(gene_ctrl_expressions) + fold_change[1]
        else:
            lower = np.array(gene_ctrl_expressions) * fold_change[0]
            upper = np.array(gene_ctrl_expressions) * fold_change[1]
        indices = np.arange(len(self.repository))
        if shuffle is True:
            np.random.shuffle(indices)
        for i in indices:
            measures = self.repository[i]
            t, p = stats.ttest_ind(lower, measures, equal_var=True)
            if t > 0 or (p/2.) > self.alpha:
                continue
            t, p = stats.ttest_ind(upper, measures, equal_var=True)
            if t < 0 or (p/2.) > self.alpha:
                continue
            return measures
        if force is False:
            msg = 'Cannot find expression measures with requested fold changes.'
            raise NonExistingExpressionException(msg)
        else:
            change = np.random.normal(loc=np.mean(fold_change), scale=std,
                                      size=len(gene_ctrl_expressions))
            if log_scale is True:
                return gene_ctrl_expressions + change
            else:
                return gene_ctrl_expressions * change



class WilcoxonRankSumDExpress(DExpress):
    """A criteria based on independent Wilcoxon rank-sum test.

    Args:
        repository (Repository): A expression repository to be used for
            simulating differential expression.
        alpha (float): A positive number between 0 and 1 used as
            significant level.

    """

    def __init__(self, repository, alpha):
        super().__init__(repository)
        self.alpha = alpha

    def __call__(self, gene_ctrl_expressions, fold_change, force=True,
                 shuffle=True, std=0.5, log_scale=True):
        """Simulate differential expression.

        Args:
            gene_ctrl_expressions(array-like): A sequence of expression values.
            fold-change (float): A fold-change value. Note that this value is
                interpreted along with the 'log_scale' parameter. Expression data are
                commonly log normalized to deal with heterogeneity of variance.
                In such cases the fold_change value should be in log scale.
            force (bool): Force the method to return an array of expression
                values. If there is no row in the repository that meet the
                criteria, gene_ctrl_expressions will be shifted using normally
                distributed values where the mean is a random value in the
                fold-change interval and a standard deviation of std.
            shuffle(bool): If True (default), rows of repository will be
                searched in random order for finding the row that meets the
                differential expression criteria; otherwise, the repository
                will be searched from the first row to the last row.
            std(float): A positive value respresenting standard deviation of
                the normally distributed shifts. Only applicable when force
                is True.
            log_scale (bool): True if expression data is in log scale, i.e. data has been
                log normalized, False otherwise.

        Returns:
            numpy.ndarray: An element of the respository that meets the
                differential expression criteria.

        Raises:
            TypeError: If fold_change is not a tuple of length two.
            ValueError: If fold_change[0] >= fold_change[1].

        """
        # Check fold_change type
        if not isinstance(fold_change, tuple) or len(fold_change) != 2:
            raise TypeError('fold_change must be a tuple of length 2.')
        if fold_change[0] > fold_change[1]:
            msg = 'fold_change[0] must be less than fold_change[1]'
            raise InvalidLogFoldChangeValue(msg)
        if fold_change[0] * fold_change[1] <= 0:
            msg = ('fold_change must be a tuple of two numbers of the '
                   'same sign.')
            raise InvalidLogFoldChangeValue(msg)
        average_fold_change = (fold_change[0] + fold_change[1]) / 2
        if log_scale is True:
            base_expression = np.array(gene_ctrl_expressions) + average_fold_change
        else:
            base_expression = np.array(gene_ctrl_expressions) * average_fold_change
        indices = np.arange(len(self.repository))
        if shuffle is True:
            np.random.shuffle(indices)
        for i in indices:
            measures = self.repository[i]
            t, p = stats.ranksums(base_expression, measures, equal_var=True)
            if p < self.alpha:
                continue
            return measures
        if force is False:
            msg = 'Cannot find expression measures with requested fold changes.'
            raise NonExistingExpressionException(msg)
        else:
            change = np.random.normal(loc=np.mean(fold_change), scale=std,
                                      size=len(gene_ctrl_expressions))
            if log_scale is True:
                return gene_ctrl_expressions + change
            else:    
                return gene_ctrl_expressions * change


class FoldChangeDExpress(DExpress):
    """A criteria based on fold-change.

    Args:
        repository (Repository): A expression repository to be used for
            simulating differential expression.
        alpha (float): A positive number between 0 and 1 used as
            significant level.

    """

    def __init__(self, repository, alpha):
        super().__init__(repository)
        self.alpha = alpha

    def __call__(self, gene_ctrl_expressions, fold_change, log_scale=True):
        """Simulate differential expression.

        Args:
            gene_ctrl_expressions(array-like): A sequence of expression values.
            fold-change (float): A fold-change value. Note that this value is
                interpreted along with the 'log_scale' parameter. Expression data are
                commonly log normalized to deal with heterogeneity of variance.
                In such cases the fold_change value should be in log scale.
            log_scale (bool): True if expression data is in log scale, i.e. data has been
                log normalized, False otherwise.


        Returns:
            numpy.ndarray: An element of the repository that meets the
                differential expression criteria.

        Raises:
            TypeError: If fold_change is not a tuple of length two.
            ValueError: If fold_change[0] >= fold_change[1].

        """
        # Check fold_change type
        if not isinstance(fold_change, tuple) or len(fold_change) != 2:
            raise TypeError('fold_change must be a tuple of length 2.')
        if fold_change[0] > fold_change[1]:
            msg = 'fold_change[0] must be less than fold_change[1]'
            raise InvalidLogFoldChangeValue(msg)
        if fold_change[0] * fold_change[1] <= 0:
            msg = ('fold_change must be a tuple of two numbers of the '
                   'same sign.')
            raise InvalidLogFoldChangeValue(msg)
        fc = np.random.uniform(low=fold_change[0], high=fold_change[1])
        if log_scale is True:
            base_expression = np.array(gene_ctrl_expressions) + fc
        else:
            base_expression = np.array(gene_ctrl_expressions) * fc
        return base_expression

