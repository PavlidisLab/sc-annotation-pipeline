# nf-core/scannotate: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v2.0.0dev - [Unreleased]

### `Added`

- Reorganized pipeline to follow nf-core best practices
- Added modular structure with local modules and subworkflows
- Added `conf/base.config` for resource management
- Added `conf/modules.config` for module-specific settings
- Added `conf/test.config` for minimal testing
- Added stub blocks to all modules for dry-run testing
- Added `versions.yml` output to all modules for software version tracking

### `Changed`

- Moved from single `sc-annotate.nf` to modular `main.nf` + `workflows/scannotate.nf`
- Renamed processes to SCREAMING_SNAKE_CASE per nf-core conventions
- Moved asset files from `meta/` to `assets/`
- Reorganized modules into `modules/local/<tool>/main.nf` structure
- Created subworkflows for logical groupings of modules

### `Fixed`

### `Dependencies`

### `Deprecated`

## v1.2.1 - 2024-XX-XX

- Previous stable release before nf-core reorganization
- Added family-level cell type annotations
- Made outlier fields configurable
- Updated QC report configuration
