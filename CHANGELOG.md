# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.5.0] - 2024-12-16
### Changed
- vdW radii scaling parameter can now be adjusted via `mindlessgen.toml` or CLI
- check_distance function now checks based on the sum of the van der Waals radii and a scaling factor acessible via `mindlessgen.toml` or CLI
- better type hints for `Callables`
- clearer differentiation between the distinct scaling factors for the van der Waals radii
- `README.md` with more detailed explanation of the element composition function
- Default `max_cycles` for the generation & refinement set to 200
- Allow fixed molecule compositions in a simpler way
- `check_config` now ConfigClass-specific
- modify atom list adaption to `element_composition` such that a random integer in the given range is taken and not the lower/upper bound

### Fixed
- unit conversion for (currenly unused) vdW radii from the original Fortran project
- minor print output issues (no new line breaks, more consistent verbosity differentiation, ...)
- bug in `postprocess_mol` which led to an unassigned return variable in the single-point case
- bug leading to `UnicodeDecodeError` when reading `xtb` output files
- bug with all atom lists being initialized with a length of 102 instead of 103
- inconsistent default values for the `mindlessgen.toml` and the `ConfigManager` class
- legacy pseudo random number generation removed and replaced by `np.random.default_rng()` for avoiding interference with other packages

### Added
- support for the novel "g-xTB" method (previous working title: GP3-xTB)
- function which contracts the coordinates after the initial generation
- function which is able to printout the xyz coordinates to the terminal similar to the `.xyz` layout
- elements 87 to 103 are accessible via the element composition. If `xtb` is the engine, the elements will be replaced by their lighter homologues.
- support for `python-3.13`
- option to set a fixed molecular charge, while ensuring `uhf = 0`
- `element_composition` and `forbidden_elements` can now be directly set to a `dict` or `list`, respectively, via API access

### Breaking Changes
- Removal of the `dist_threshold` flag and in the `-toml` file.
- The number of unpaired electrons (`Molecule.uhf`) is now set to 0 if `xtb` is used as `QMMethod` and a lanthanide is within the molecule to match the `f-in-core` approximation.
- "Contract Coordinates" functionality set to `true` by default in the `mindlessgen.toml` file.
- `basename.UHF` and `basename.CHRG` are only written to disk if they differ from the default value (0 and 0, respectively).

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
