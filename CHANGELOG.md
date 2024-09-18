# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]
### Changed
- Default file name of `.xyz` file contains prefix `mlm_`
- Comment line of `.xyz` file contains the total charge and number of unpaired electrons
- Default ORCA calculation changed from r2SCAN-3c to PBE/def2-SVP
- `verbosity = 3` always prints full QM output
- Adapted generation of number of unpaired electrons; thereby, support for Ln's
- Shifted group / element sorting definitions to miscellaneous
- `xyz` files are written on the fly, and not post-generation
- `forbidden_elements` and `element_composition` influences hydrogen and organic element addition
- GFN<n>-xTB level can now be set

### Added
- Optimization via DFT in the post-processing step
- Detailed input of ORCA settings (functional, basis, grid size, SCF cycles, ...) possible
- Maximum number of optimization cycles are an argument for the `QMMethod.optimize` base function
- Debug option for the refinement and post-processing step specifically
- Return type for `single_molecule_generator`
- Check for consistency of the `min_num_atoms` and `max_num_atoms` constraint
- Similar to the `<basename>.CHRG` file, also a `<basename>.UHF` is printed
- HOMO-LUMO gap check within the refinement step and corresponding Config option called "refine_hlgap"
- `GeneralConfig` switch for writing `xyz` files

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
