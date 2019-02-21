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

## Versioning

We use [Semantic Versioning 2.0.0](http://semver.org/) for versioning.

## Authors

* [**Farhad Maleki**](https://github.com/FarhadMaleki) - *Initial work* 

