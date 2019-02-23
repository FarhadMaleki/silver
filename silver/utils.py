"""This module provides utility functions for silver.

"""

import numpy as np
import pandas as pd
from silver.dataset import Dataset
from silver.exceptions import InvalidContrastException
from silver.expression_profile import ExpressionProfile
from silver.repository import Repository
from silver.dexpress import TTestDExpress


def read_contrast(address, ctrl_symbol='c', case_symbol='d', sep='\t'):
    """Read contrast from a file.

    Args:
        address (str): Address of a file that determines the phenotype
            classes.
        ctrl_symbol (str): Symbol used in contrast file for representing
            control samples. Default value is 'c'.
        case_symbol (str): Symbol used in contrast file for representing
            case samples. Default value is 'd'.
        sep (str): Field separator.  Default value is '\t'.

    Returns:
        tuple: A tuple of size two, which contains control indices and
            case indices, respectively.
    """
    contrast = []
    with open(address, 'r') as fin:
        for line in fin:
            line = line.strip()
            if line == "":
                continue
            contrast.extend([phenotype.strip()
                             for phenotype in line.split(sep)])
    if ctrl_symbol not in contrast:
        raise InvalidContrastException('contrast must contain ctrl_symbol' +
                                       ' ({}).'.format(ctrl_symbol))
    num_samples = len(contrast)
    contrast = np.array(contrast)
    ctrl_col_indices = np.arange(num_samples)[contrast == ctrl_symbol]
    case_col_indices = np.arange(num_samples)[contrast == case_symbol]

    return list(ctrl_col_indices), list(case_col_indices), list(contrast)


###############################################################################
def read_profile(address, **kwargs):
    """ Read an expression profile.

    Args:
        address (str): Address of expression profile.
        **kwargs: The as kwargs in pandas.read_csv.

    Returns:
        ExpressionProfile: The expression profile read from the address.

    """
    df = pd.read_csv(address, **kwargs)
    df.index = [str(gene_id) for gene_id in df.index]
    return ExpressionProfile(df)


###############################################################################
def read_fold_change_file(address, id_col_name, lower_bound_col,
                          upper_bound_col, sep='\t'):
    """Read fold change file to a dictionary.

    Args:
        address (str): The address of the fold change file, which is a
            tabular file with three named columns.
        id_col_name (str): The name of the column that contains gene
            names/ids.
        lower_bound_col (str): The name of the column that contains
            lower bounds for fold changes.
        upper_bound_col (str): The name of the column that contains
            upper bounds for fold changes.
        sep (str): The field separator in fold change file.
            The default is '\t'.

       Returns:
        dict: A dictionary where gene names/ids are keys and tuples of
            (lower bound, upper bound) for gene fold changes as values.

    """
    df = pd.read_csv(address, index_col=False, sep=sep)
    df = df[[id_col_name, lower_bound_col, upper_bound_col]]
    df[id_col_name] = [str(e) for e in df[id_col_name]]
    gfc = {}
    for _, row in df.iterrows():
        gfc[row[id_col_name]] = (row[lower_bound_col], row[upper_bound_col])
    return gfc


###############################################################################
def simulate(profile_address=None,
             contrast_address=None,
             gene_fc_address=None,
             fc_id_col_name=None,
             lower_bound_col_name=None,
             upper_bound_col_name=None,
             num_simulated_ctrls=None,
             num_simulated_cases=None,
             profile_id_col_name=None,
             num_repository_reps=1,
             ctrl_symbol='c',
             case_symbol='d',
             profile_sep='\t',
             contrast_sep='\t',
             fold_change_sep='\t',
             alpha=0.05,
             random_state=0):
    """

    Args:
        profile_address (str): Address of a real expression profile.
            This file can be downloaded from online sources such as
            GEO and ArrayExpress.
        contrast_address (str): Address of a file that determine
            the phenotype of each sample in the real expression.
        gene_fc_address (str): Address of a tabular file with three
            named columns. This file contains the genes that a user
            wishes to differentially express. Each row of this file
            contains the name/ID of a gene from the expression profile,
            and the lower and upper bound for fold change between case
            and control samples. Down-regulation must be represented
            with negative numbers and up-regulation with positive numbers.
        fc_id_col_name (str): The column name for the gene name/Id in
            fold change file.
        lower_bound_col_name (str): The column name for the lowest fold
            change allowed for each gene.
        upper_bound_col_name (str): The column name for the highest fold
            change allowed for each gene.
        num_simulated_ctrls (int): Number of control samples to be
            simulated. It must be less than or equal to the number of
            original control samples.
            Also, num_simulated_ctrls + num_simulated_cases must be less
            than or equal to the number of actual controls samples.
        num_simulated_cases (int): Number of case samples to be
            simulated. It must be less than or equal to the number of
            original case samples.
            Also, num_simulated_ctrls + num_simulated_cases must be less
            than or equal to the number of actual controls samples.
        profile_id_col_name (str): The column name for the gene name/Id in
            expression profile.
        num_repository_reps (int): The number of subsets of
            size num_simulated_cases that are used during differential
            expression. The default is 1.
        ctrl_symbol (str): The string used to represent a control sample
            in contrast file. The default is 'c'.
        case_symbol (str): The string used to represent a case sample
            in contrast file. The default is 'd'.
        profile_sep (str): Field separator.  Default value is '\t'.
        contrast_sep (str): Field separator.  Default value is '\t'.
        fold_change_sep (str): Field separator.  Default value is '\t'.
        alpha (float): Significance level. The default is 0.05.
        random_state (int): The default is 0.
    """
    #Set random seed
    np.random.seed(random_state)
    #Make sure that information about input file has been entered
    msg = 'Information about input files must be provided.'
    assert None not in {profile_address, contrast_address, gene_fc_address,
                        fc_id_col_name, lower_bound_col_name,
                        upper_bound_col_name, num_simulated_ctrls,
                        num_simulated_cases, profile_id_col_name}, msg
    # Read contrast from a file
    temp = read_contrast(contrast_address, ctrl_symbol=ctrl_symbol,
                         case_symbol=case_symbol, sep=contrast_sep)
    ctrl_indices, case_indices, contrast = temp
    # Read the original expression profile from a file
    profile = read_profile(profile_address,
                           sep=profile_sep,
                           index_col=profile_id_col_name)
    # Read fold change information from a file
    gfc = read_fold_change_file(gene_fc_address, fc_id_col_name,
                                lower_bound_col_name, upper_bound_col_name,
                                sep=fold_change_sep)
    # Instantiate a Dataset object
    dataset = Dataset(profile.samples(ctrl_indices),
                      profile.samples(case_indices))
    # Generate simulated controls and a replicates (no differential expression)
    sim_ctrls, sim_cases = dataset.make_replicates(num_simulated_ctrls,
                                                   num_simulated_cases)
    # Assemble a repository
    repository = Repository(profile.samples(case_indices), num_simulated_cases,
                            num_repetitions=num_repository_reps,
                            random_state=random_state)
    # Create a DEpress object
    dexpress_obj = TTestDExpress(repository, alpha=alpha)
    # Apply differential expression
    dataset.diff_express(sim_ctrls, sim_cases, gfc, dexpress_obj)
    return sim_ctrls, sim_cases
