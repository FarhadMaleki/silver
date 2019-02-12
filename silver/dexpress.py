"""This module contains criteria used for simulating differential expression.

Each criterion must extend DExpress class, which is an abstract class.

"""
from abc import ABC, abstractmethod


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
