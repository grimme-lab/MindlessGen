# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]
### Changed
- Default file name of `.xyz` file contains prefix `mlm_`
- Comment line of `.xyz` file contains the total charge and number of unpaired electrons
- Default ORCA calculation changed from r2SCAN-3c to PBE/def2-SVP

### Added
- Optimization via DFT in the post-processing step
- Detailed input of ORCA settings (functional, basis, grid size, SCF cycles, ...) possible
- Maximum number of optimization cycles are an argument for the `QMMethod.optimize` base function

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
