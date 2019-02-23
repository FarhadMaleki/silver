# Silver

Silver is a framework for the evaluation of gene set analysis methods. It provides a mechanism for simulating gene expression datasets while avoiding oversimplifying assumptions about the distribution of gene expression values. Gene expression data is simulated using real expression datasets (as input) while having control of differential expression of genes. Silver also makes a quantitative evaluation of gene set analysis methods possible.

## Getting Started
The following steps will guide you through the process of running Silver on your local machine.

### Prerequisites

The main dependencies of Silver are numpy, pandas, and scipy. See _**requirements.txt**_ for the list of requirements.

### Installing
It is a good practice to use a virtual environment for deploying Python programs. Using **conda**, we will create an environment named *Silver*. The environment name is arbitrary.

```bash
conda create -n silver python=3.6
```
To install requirements, the following command can be run.

```bash
make setup
```

## Running the tests

Regression tests can be run through the following command:

```bash
make regression
```

### Adding more tests

New tests should be added as modules where their names start with *test_* under *test* directory.

## Documentation

A comprehensive documentation is available [here](https://farhadmaleki.github.io/silver/).

## Example

For convenience, the following example is provided as a python script named **usecase.py**.  You should be able to run the code in Python 3 as follows:
``` bash
python usecase.py
```
Note: Make sure that **usecase.py** is in your current directory.

---

To simulate a dataset, we first need to import the following modules from ```silver```.

``` python
from silver.utils import simulate
```

The follwoig files are required to simulate a dataset.
<dl>
    <dt>*A real expression profile*</dt>
    <dd>This file can be downloaded from online sources such as GEO and ArrayExpress.</dd>
    <dt>*A contrast file*</dt>
    <dd>This detrmines the phenotype of each sample in the real expression.</dd>
    <dt>*A fold change file*</dt>
    <dd> This file contains the genes that a user wishes to differentially 
express. Each row of this file contains the name/ID of a gene from the expression profile, and the lower and upper bound for fold change between case and control samples. Down-regulation must be represented with negative numbers and up-regulation with positive numbers.</dd>
</dl>

``` python
PROFILE_ADDRESS = "data/GSE53757_profile.txt"
CONTRAST_ADDRESS = "data/GSE53757_contrast.txt"
GENE_FC_ADDRESS = "data/DE_gene_fold_change.txt"
```

See the *data* directory for an example for each of these files. 
Then we need to read these files. To do so, we need the following information about each of them:

<dl>
    <dt>*FC_ID_COL_NAME*</dt>
    <dd>The column name for the gene name/Id in fold change file.</dd>
    <dt>*LOWER_BOUND_COL_NAME*</dt>
    <dd>The column name for the lowest fold change allowed for each gene.</dd>
    <dt>*UPPER_BOUND_COL_NAME*</dt>
    <dd>The column name for the highest fold change allowed for each gene.</dd>
    <dt>PROFILE_ID_COL_NAME</dt>
    <dd>The column name for the gene name/Id in expression profile.</dd>
</dl>

These are as follows for our example:

``` python
LOWER_BOUND_COL_NAME = 'FCLower'
UPPER_BOUND_COL_NAME = 'FCUpper'
PROFILE_ID_COL_NAME = 'ID'
FC_ID_COL_NAME = 'ID'
```

Then we set the number of simulated controls and simulated cases. 

``` python
NUM_SIMULATED_CTRLS = 20
NUM_SIMULATED_CASES = 20
```

Using this information, we are now able to simulate a dataset using simulate_data
from Silver.

``` python 
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
```

Finally, we will combine the simulated controls and cases, get the expression profile and write the results to a file.

``` python
expression_dataset = sim_ctrls.concat(sim_cases)
result = expression_dataset.profile
simulated_profile_address = 'simulated.profile.txt'
result.to_csv(simulated_profile_address, sep='\t')
```

We can also create a contrast file for the simulated dataset as follows:

``` python 
simulated_contrast_address = 'simulated.contrast.txt'
with open(simulated_contrast_address, 'w') as fout:
    fout.write('{}\n'.format('\t'.join(['c'] * NUM_SIMULATED_CTRLS +
                                       ['d'] * NUM_SIMULATED_CASES)))
```

## Versioning

We use [Semantic Versioning 2.0.0](http://semver.org/) for versioning.

## Authors

* [**Farhad Maleki**](https://github.com/FarhadMaleki) - *Initial work* 

