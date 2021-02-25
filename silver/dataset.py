"""This module simulates expression datasets.

"""
import random
import logging
import numpy as np
from silver.exceptions import MismatchedProfileException
from silver.exceptions import NonExistingExpressionException


logger = logging.getLogger(__name__)

class Dataset(object):
    """ Set the required case/control data to be used for Dataset creation.
    Args:
        ctrl_profile (ExpressionProfile): A ExpressionProfile object
            containing expression measures for control samples.
        case_profile (ExpressionProfile): A ExpressionProfile object
            containing expression measures for case samples.

    Raises:
        MismatchedProfileException: If the order of gene ids in
            ctrl_profile and case_profile are not the same.

    """
    def __init__(self, ctrl_profile, case_profile):
        if not np.array_equal(list(ctrl_profile.keys()),
                              list(case_profile.keys())):
            msg = ('The order of gene ids in ctrl_profile and case_profile ' +
                   'must be the same.')
            raise MismatchedProfileException(msg)
        self.controls = ctrl_profile
        self.cases = case_profile

    ############################################################################
    def make_custom_replicates(self, simulated_ctrl_indices,
                               simulated_case_indices,
                               replace=False):
        """Create a replicate datasets with predetermined samples.

        Args:
            simulated_ctrl_indices (array-like): The indices of samples
                to be used as control samples. These samples are chosen from
                all controls used for initializing the object.
            simulated_case_indices (array-like): The indices of samples
                to be used as case samples. These samples are chosen from
                all controls used for initializing the object.
            replace (bool): True for using sampling with replacement and False
                otherwise. Default is False.

        Raises:
            IndexError: If simulated_ctrl_indices (simulated_case_indices)
                contains a value greater than or equal to the total number
                of controls (cases).
            ValueError: If elements of simulated_ctrl_indices
                or simulated_case_indices are not unique, or if there is an
                element that is in both simulated_ctrl_indices and
                simulated_case_indices. Only applicable when the replace
                parameter is False.

        Returns:
            (tuple): simulated_ctrls and simulated_cases, where simulated_cases
                are replicate (with no differential expression) for
                simulated_ctrls.

        """
        if replace is False:
            msg = '{} must contain unique values.'
            if len(set(simulated_ctrl_indices)) != len(simulated_ctrl_indices):
                raise IndexError(msg.format('simulated_ctrl_indices'))
            if len(set(simulated_case_indices)) != len(simulated_case_indices):
                raise IndexError(msg.format('simulated_case_indices'))
            if (set(simulated_ctrl_indices) & set(simulated_case_indices)) != set():
                raise ValueError('Group indices must not overlap.')

        simulated_ctrls = self.controls.samples(simulated_ctrl_indices)
        simulated_cases = self.controls.samples(simulated_case_indices)
        return simulated_ctrls, simulated_cases

    ###########################################################################
    def make_replicates(self, num_sim_ctrls, num_sim_cases,
                         replace=False, random_state=None):
        """Create a replicate dataset with predetermined number of samples.

        Args:
            num_sim_ctrls (int): Number of control samples to be simulated.
            num_sim_cases (int): Number of case samples to be simulated.
            replace (bool): True for using sampling with replacement and False
                otherwise. Default is False.
            random_state (int): A seed for reproducing method results.

        Raises:
            IndexError: If simulated_ctrl_indices (simulated_case_indices)
                contains a value greater than or equal to the total number
                of controls (cases).
            ValueError: If elements of simulated_ctrl_indices
                or simulated_case_indices are not unique, or if there is an
                element that is in both simulated_ctrl_indices and
                simulated_case_indices.

        Returns:
            (tuple): simulated_ctrls and simulated_cases, where simulated_cases
                are replicate (with no differential expression) for
                simulated_ctrls.
        """
        if random_state is not None:
            np.random.seed(random_state)
        total_num_controls = self.controls.shape[1]
        if replace is False and (num_sim_ctrls + num_sim_cases > total_num_controls):
            raise ValueError('num_sim_ctrls + num_sim_cases must ' +
                             'be less than or equal to the total number of ' +
                             'controls. Alternatively you should set the "replace" ' +
                             'parameter as True to use sampling with replacement ' +
                             'instead of sampling without replacement.')
        indices = np.arange(total_num_controls)
        if replace is True:
            indices = random.choices(indices, k=num_sim_ctrls + num_sim_cases)
        indices = np.random.permutation(indices)
        simulated_ctrl_indices = indices[:num_sim_ctrls]
        simulated_case_indices = indices[num_sim_ctrls:(num_sim_ctrls + num_sim_cases)]
        return self.make_custom_replicates(simulated_ctrl_indices,
                                           simulated_case_indices,
                                           replace=replace)

    ############################################################################
    @staticmethod
    def diff_express(ctrls, cases, genes_fold_changes, dexpress_obj):
        """Differentially express some genes on the simulated cases.

        Args:
            ctrls (ExpressionProfile): A ExpressionProfile containing simulated
                controls.
            cases (ExpressionProfile): A ExpressionProfile containing simulated
                cases, which are
                a replicate for controls without differential expression.
            genes_fold_changes (dict): A dictionary representing the
                fold-change values for some genes. Each gene in
                genes_fold_changes.keys() must be an index element of
                ctrls and cases.
            dexpress_obj (DExpress): A DExpress object that determines the
                procedure to be applied for differential expression.

        Returns:
            A tuple containing ctrls (unchanged) and updated cases.

        """
        for gene, fold_change in genes_fold_changes.items():
            expressions = ctrls.get(gene)
            try:
                expressions = dexpress_obj(expressions, fold_change)
            except NonExistingExpressionException:
                msg = ('{} cannot be expressed by' +
                       ' {} fold change').format(gene, fold_change)
                logger.warning(msg)
            finally:
                cases.set(gene, expressions)
        return ctrls, cases

