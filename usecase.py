"""This module shows an example for simulating expression profiles.

"""
from silver.utils import simulate


# Input informaiton
PROFILE_ADDRESS = "data/GSE53757_profile.txt"
CONTRAST_ADDRESS = "data/GSE53757_contrast.txt"
GENE_FC_ADDRESS = "data/DE_gene_fold_change.txt"
# Input column Names
LOWER_BOUND_COL_NAME = 'FCLower'
UPPER_BOUND_COL_NAME = 'FCUpper'
PROFILE_ID_COL_NAME = 'ID'
FC_ID_COL_NAME = 'ID'
# Set the number of simulated controls and simulated cases
NUM_SIMULATED_CTRLS = 20
NUM_SIMULATED_CASES = 20
# Simulate control and case samples
sim_ctrls, sim_cases = simulate(profile_address=PROFILE_ADDRESS,
                                contrast_address=CONTRAST_ADDRESS,
                                gene_fc_address=GENE_FC_ADDRESS,
                                fc_id_col_name=FC_ID_COL_NAME,
                                lower_bound_col_name=LOWER_BOUND_COL_NAME,
                                upper_bound_col_name=UPPER_BOUND_COL_NAME,
                                num_simulated_ctrls=NUM_SIMULATED_CTRLS,
                                num_simulated_cases=NUM_SIMULATED_CASES,
                                profile_id_col_name=PROFILE_ID_COL_NAME,
                                num_repository_reps=10,
                                ctrl_symbol='c',
                                case_symbol='d',
                                profile_sep='\t',
                                contrast_sep='\t',
                                fold_change_sep='\t',
                                alpha=0.05,
                                random_state=123456)
# Combine the simulated controls and cases
expression_dataset = sim_ctrls.concat(sim_cases)
# Get the expression profile
result = expression_dataset.profile
# Write the simulated dataset to a file
simulated_profile_address = 'simulated.profile.txt'
result.to_csv(simulated_profile_address, sep='\t')
# Make a contrast for the simulated dataset and write it to a file
simulated_contrast_address = 'simulated.contrast.txt'
with open(simulated_contrast_address, 'w') as fout:
    fout.write('{}\n'.format('\t'.join(['c'] * NUM_SIMULATED_CTRLS +
                                       ['d'] * NUM_SIMULATED_CASES)))
