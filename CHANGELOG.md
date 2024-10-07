# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]
### Changed
- vdW radii scaling parameter can now be adjusted via `mindlessgen.toml` or CLI
- The check_distance function now checks based on the sum of the van der Waals radii and a scaling factor acessible via `mindlessgen.toml` or CLI
- better type hints for `Callables`

### Fixed
- Unit conversion for (currenly unused) vdW radii from the original Fortran project
- minor print output issues (no new line breaks, more consistent verbosity differentiation, ...)
- bug in `postprocess_mol` which led to an unassigned return variable in the single-point case

### Added
- Support for the novel "g-xTB" method (working title: GP3-xTB)
- A function which contracts the coordinates after the initial generation to enhance the convergence with the qm method.

### Breaking Changes
- Removal of the `dist_threshold` flag and in the `-toml` file.

### Breaking Changes
- Removal of the `dist_threshold` flag and in the `-toml` file.

## [0.4.0] - 2024-09-19
### Changed
- Default file name of `.xyz` file contains prefix `mlm_`
- Comment line of `.xyz` file contains the total charge and number of unpaired electrons
- Default ORCA calculation changed from r2SCAN-3c to PBE/def2-SVP
- `verbosity = 3` always prints full QM output
- Adapted generation of number of unpaired electrons; thereby, support for Ln's
- Shifted group / element sorting definitions to miscellaneous
- `xyz` files are written on the fly, and not post-generation
- GFN<n>-xTB level can now be set
- `mindless.molecules` file is written continuously during generation

### Fixed
- `test_iterative_optimization` more deterministic
- wrong atom range check in for the isomerization mode ([#21](https://github.com/grimme-lab/MindlessGen/pull/21))
- `forbidden_elements` and `element_composition` influences hydrogen and organic element addition
- more realistic default `mindlessgen.toml` entries

### Added
- Optimization via DFT in the post-processing step
- Detailed input of ORCA settings (functional, basis, grid size, SCF cycles, ...) possible
- `min_num_atoms` and `max_num_atoms` consistency check
- Maximum number of optimization cycles are an argument for the `QMMethod.optimize` base function
- Debug option for the refinement and post-processing step specifically
- Return type for `single_molecule_generator`
- Check for consistency of the `min_num_atoms` and `max_num_atoms` constraint
- Similar to the `<basename>.CHRG` file, also a `<basename>.UHF` is printed
- HOMO-LUMO gap check within the refinement step and corresponding Config option called "refine_hlgap"
- `GeneralConfig` switch for writing `xyz` files
- `PyPi` and `TestPyPi` upload of releases (new workflow)

## [0.3.0] - 2024-08-20
### Breaking Changes
- ...

### Added
- ...

### Changed
- ...

### Removed
- ...

### Fixed
- ...

### Deprecations
- ...
