# nf-core/scannotate: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v2.1.1 - 2026-10-19

### `Fixed`

- The mouse reference lost two of its three Census datasets (and the Macrophage, Pericyte, Neural stem cell, OPC and CA2-IG-FC subclasses) when their gitignored author label files were absent; it now matches the benchmark reference (26 subclasses).
- The reference is subsampled per renamed subclass rather than per Census `cell_type`, as in the benchmark (applies to human and mouse).

### `Added`

- `assets/REFERENCE_INPUTS.sha256` and `scripts/reference_inputs.sh` record and check the SHA-256 of the reference build inputs, including gitignored files.

### `Changed`

- Removed author label files and `original_celltype_columns.tsv` rows for datasets outside the default `ref_collections`.

## v2.1.0 - 2026-10-15

### `Added`

- Built the CTA protocol from the pipeline configuration: `sc-pipeline-<version>-<classifier>-census<census>-ref<tag>-n<subsample_ref>-cutoff<cutoff>`, with `cta_protocol` (override) and `ref_tag` parameters and `assets/ref_tags.tsv`
- Created the CTA protocol in Gemma from the pipeline (`addProtocol` in `LOAD_CTA`), with a description recording the configuration; an existing protocol is not an error.
- Added `--classifier knn` option for the scVI embedding classifier (`rf` stays the default); process `RF_CLASSIFY` renamed to `CLASSIFY_SCVI`
- Added `use_gemma` master switch, with `upload_gemma` nested under it
- Added full human and mouse test profile matrix and `scripts/run_test_profiles.sh`
- Added support for a pre-combined h5ad as `study_path`
- Added Ma et al. reference collection, rename mappings and markers for human immune, vascular and erythroid cell types
- Added class and family marker heatmaps and QC outlier mask and tables to MultiQC
- Added Apache-2.0 license, `audit_docs/`, and a script to clean stale-version Gemma CTAs (`scripts/cleanup_stale_ctas.py`)

### `Changed`

- Replaced the workflow diagram with a Mermaid flow diagram (`docs/diagrams/scannotate_flow.mmd`)
- Dropped Tabula Muris Senis from mouse `ref_collections` and mouse `subsample_ref` 100
- Mouse rename map gained `OPC` and `L6 IT Car3` rows
- Merged two human VIP subclass names into `VIP GABAergic interneuron` and remapped `Micro-PVM` to `Microglia`
- `scvi_classify.py` now requires `--ref_path`
- Added container and scVI model caching, and SLURM CPU tiers by process label
- Split `bin/utils.py` into classify, QC, marker and census submodules
- Removed hardcoded Gemma credentials from script defaults
- Gemma credentials are read from the environment (`GEMMA_USERNAME`, `GEMMA_PASSWORD`) and no longer appear on command lines or in work files
- Shortened the README to the essentials; the full text is in git history

### `Fixed`

- `use_staging` now honored in the Gemma metadata fetch and the production gemmapy path
- Fixed sample metadata misalignment in `get_gemma_meta.py`
- Fixed apptainer container resolution
- Fixed MEX ingestion: categorical `sample_name` crash, unrecoverable samples now skipped, blank Ensembl ids dropped
- Fixed dead paths and broken test fixtures
- Cells with zero counts are removed before QC and classification (a sample that loaded empty crashed the QC fit); the number removed per sample is logged
- `get_gemma_meta.py` falls back to the dataset taxon when the platform has none (generic platforms such as `Generic_mouse_ncbiIds` returned a null taxon)

## v2.0.0 - 2026-01-12

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
