# Mindless Molecule Generator

<img style="float: right;" src="assets/C1H2N1O2Te2Er1Lu2_89bd3e.png" width="300">

![CI](https://github.com/marcelmbn/MindlessGen/actions/workflows/ci.yml/badge.svg)
<a href="http://www.apache.org/licenses/LICENSE-2.0">
  <img src="https://img.shields.io/badge/License-Apache%202.0-orange.svg" alt="Apache-2.0"/>
</a>
<a href="https://img.shields.io/badge/Python-3.10%20|%203.11%20|%203.12-blue.svg">
  <img src="https://img.shields.io/badge/Python-3.10%20|%203.11|%203.12-blue.svg" alt="Python Versions"/>
</a>

## Installation

> [!IMPORTANT]  
> `xtb` has to be available on your machine, either via a `conda-forge` installation, a release binary, or compiled from source. Further information is available [here](https://github.com/grimme-lab/xtb).

### Non-development purposes

The project can be simply installed in an existing virtual environment (by, e.g., `conda` or `mamba` (see also [here](https://github.com/conda-forge/miniforge) and [here](https://conda.io/projects/conda/en/latest/user-guide/getting-started.html))) with
```
pip install .
```

A matching Python environment can be set up and activated via the following command using the tools above:
```
mamba create -n mindlessgen python=3.12
mamba activate mindlessgen
```

### Development purposes

For working on the code of `mindlessgen`, the following setup is recommended:
```
mamba create -n mindlessgen python=3.12
mamba activate mindlessgen
pip install -e '.[dev]'
```
Thereby, all necessary development tools (e.g., `ruff`, `mypy`, and `pre-commit`) are installed.
Before start to make changes in the code, activate the `pre-commit` hooks via:
```
pre-commit install
```
Before pushing a commit to the repository, please run also the optional tests, which depend on external dependencies like xtb, via
```
pytest -vv --optional
```

## Usage

> [!WARNING]
> `mindlessgen` is still subject of drastic API changes and in the early development phase.

`mindlessgen` can be executed after installation in the desired environment via:
```
mindlessgen -h
```
All relevent command-line options are displayed in the terminal.

When using the program for academic purposes, please cite:

_J. Chem. Theory Comput._ 2009, **5**, 4, 993–1003
```
@article{doi:10.1021/ct800511q,
author = {Korth, Martin and Grimme, Stefan},
title = {“Mindless” DFT Benchmarking},
journal = {Journal of Chemical Theory and Computation},
volume = {5},
number = {4},
pages = {993-1003},
year = {2009},
doi = {10.1021/ct800511q},
note ={PMID: 26609608},
URL = {https://doi.org/10.1021/ct800511q},
eprint = {https://doi.org/10.1021/ct800511q}
}
```


## Acknowdledgements

[T. Gasevic](https://github.com/gasevic) for creating an initial `GitHub` [migration](https://github.com/gasevic/mlmgen) of the code and providing important adjustments to the workflow.
[S. Grimme](https://www.chemie.uni-bonn.de/grimme/de/grimme) and M. Korth for the original code written in Fortran associated to the publication in [J. Chem. Theory Comput.](https://pubs.acs.org/doi/full/10.1021/ct800511q).
