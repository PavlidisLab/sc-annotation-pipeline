# nf-core/scannotate: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v2.0.0dev - [Unreleased]

### `Added`

- Added samplesheet input support (`--input`) following nf-core conventions
  - CSV format with columns: `sample`, `study_name`, `study_path`
  - Allows mixing studies to download and local paths in a single run
  - Added `assets/schema_input.json` for samplesheet validation
  - Added example samplesheet at `assets/samplesheet.csv`
- Reorganized pipeline to follow nf-core best practices
- Added modular structure with local modules and subworkflows:
  - `INPUT_CHECK` - validates inputs and downloads/prepares studies
  - `PREPARE_REFERENCE` - prepares scVI model and Census reference data
  - `PROCESS_QUERIES` - processes query data through scVI model
  - `CLASSIFY_CELLTYPES` - random forest classification
  - `QC_REPORTING` - QC analysis and MultiQC report generation
  - `GEMMA_UPLOAD` - uploads results to Gemma (optional)
- Added `conf/base.config` for resource management
- Added `conf/modules.config` for module-specific settings
- Added `conf/test.config` for minimal testing
- Added stub blocks to all modules for dry-run testing
- Added `versions.yml` output to all modules for software version tracking
- Added `ref_keys` parameter for specifying reference keys in classification
- Added `preferredCtaLevel` parameter for preferred cell type annotation level
- Added granular Gemma upload controls: `upload_cta`, `upload_clc`, `upload_multiqc`

### `Changed`

- Moved from single `sc-annotate.nf` to modular `main.nf` + `workflows/scannotate.nf`
- Renamed processes to SCREAMING_SNAKE_CASE per nf-core conventions
- Moved asset files from `meta/` to `assets/`
- Reorganized modules into `modules/local/<tool>/main.nf` structure
- Created subworkflows for logical groupings of modules
- GEMMA_UPLOAD subworkflow now has granular control over which outputs to upload

### `Fixed`

### `Dependencies`

### `Deprecated`

- `--study_names` and `--study_paths` parameters are now legacy; use `--input` samplesheet instead

## v1.2.1 - 2024-XX-XX

- Previous stable release before nf-core reorganization
- Added family-level cell type annotations
- Made outlier fields configurable
- Updated QC report configuration
