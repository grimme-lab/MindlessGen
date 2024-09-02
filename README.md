# Mindless Molecule Generator

![CI](https://github.com/grimme-lab/MindlessGen/actions/workflows/ci.yml/badge.svg)
<a href="http://www.apache.org/licenses/LICENSE-2.0">
  <img src="https://img.shields.io/badge/License-Apache%202.0-orange.svg" alt="Apache-2.0"/>
</a>
<a href="https://img.shields.io/badge/Python-3.10%20|%203.11%20|%203.12-blue.svg">
  <img src="https://img.shields.io/badge/Python-3.10%20|%203.11|%203.12-blue.svg" alt="Python Versions"/>
</a>
<img align="right" src="assets/C1H2N1O2Te2Er1Lu2_89bd3e.png" height="150" />

`mindlessgen` is a Python-based program for semi-automated generation of "mindless" small molecules, as described [here](https://pubs.acs.org/doi/full/10.1021/ct800511q).
The rule-based algorithm places atoms randomly within the coordinate space and applies several optimization, fragment detection, and sanity check steps. The program is mainly controlled via a [TOML](https://github.com/grimme-lab/MindlessGen/blob/main/mindlessgen.toml) configuration file, see below for details.

## Installation

> [!IMPORTANT]
> `xtb` (see [here](https://github.com/grimme-lab/xtb)) has to be available on your machine, either via a `conda-forge` installation, a release binary, or compiled from source. If post-processing with DFT is desired, also `orca` (see [here](https://www.faccts.de/docs/orca/6.0/manual/index.html)) has to be available.

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
Thereby, all necessary development tools (e.g., `ruff`, `mypy`, `tox`, `pytest`, and `pre-commit`) are installed.
Before start to make changes in the code, activate the `pre-commit` hooks via:
```
pre-commit install
```
Before pushing a commit to the repository, please run also the optional tests, which depend on external dependencies like `xtb`, via
```
pytest -vv --optional
```
Further information on how to contribute to this project can also be found in the [contribution guidelines](https://github.com/grimme-lab/MindlessGen/blob/main/CONTRIBUTING.md).

## Usage

> [!WARNING]
> `mindlessgen` may still be subject to API changes.

`mindlessgen` can be executed after installation in the desired environment via:
```
mindlessgen -h
```
This command displays all command line options in the terminal.
In addition, all commands are accessible via the [TOML](https://github.com/grimme-lab/MindlessGen/blob/main/mindlessgen.toml) configuration file.
The template configuration file in the root directory of the repository contains comprehensive explanations for each of the available configuration keys.
If its path is not given via `-c/--config`, the configuration file `"mindlessgen.toml"` is searched in the following paths in ascending order:
1. Current working directory (e.g., `$CWD`)
2. Home directory (e.g., `$USER/`)

The active configuration can be printed using `--print-config`.

When using the program for academic purposes, please cite:

_J. Chem. Theory Comput._ 2009, **5**, 4, 993–1003

or in `BibTeX` format:
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

## One-page overview

![One-pager overview](assets/MindlessGen_SingleSlide.png)

## Acknowdledgements

[T. Gasevic](https://github.com/gasevic) for creating an initial `GitHub` [migration](https://github.com/gasevic/mlmgen) of the code and providing important adjustments to the workflow.
[S. Grimme](https://www.chemie.uni-bonn.de/grimme/de/grimme) and M. Korth for the original code written in Fortran associated to the publication in [J. Chem. Theory Comput.](https://pubs.acs.org/doi/full/10.1021/ct800511q).
