# nf-core Reorganization Plan for sc-annotation-pipeline

## Overview

This plan outlines the changes needed to reorganize the pipeline according to [nf-core best practices](https://nf-co.re/docs/tutorials/adding_a_pipeline/overview) for eventual publication.

---

## 1. Directory Structure Changes

### Current Structure
```
cell_annotation_cortex.nf/
├── sc-annotate.nf              # Main workflow (mixed processes + workflow)
├── nextflow.config
├── params.mm.json / params.hs.json
├── bin/
├── meta/
├── modules/
│   ├── processes/
│   └── subworkflows/
├── tests/
└── README.md
```

### Target nf-core Structure
```
nf-core-scannotate/
├── main.nf                      # Entry point (minimal, calls workflows/)
├── nextflow.config              # Reorganized with includes
├── nextflow_schema.json         # NEW: Parameter schema (auto-generated)
├── modules.json                 # NEW: Track nf-core module versions
├── .nf-core.yml                 # NEW: nf-core tools configuration
│
├── conf/
│   ├── base.config              # CPU/memory requirements
│   ├── modules.config           # Module-specific settings
│   ├── test.config              # Minimal test configuration
│   ├── test_full.config         # Full test configuration
│   
│
├── workflows/
│   └── scannotate.nf            # Main pipeline workflow logic
│
├── subworkflows/
│   ├── local/
│   │   ├── download_studies/
│   │   │   └── main.nf
│   │   ├── process_queries/
│   │   │   └── main.nf
│   │   ├── classify_celltypes/
│   │   │   └── main.nf
│   │   └── qc_reporting/
│   │       └── main.nf
│   └── nf-core/                 # (for any shared nf-core subworkflows)
│
├── modules/
│   ├── local/
│   │   ├── setup_scvi/
│   │   │   ├── main.nf
│   │   │   ├── meta.yml
│   │   │   └── environment.yml
│   │   ├── get_census_adata/
│   │   │   ├── main.nf
│   │   │   ├── meta.yml
│   │   │   └── environment.yml
│   │   ├── process_query_sample/
│   │   │   └── ...
│   │   ├── process_query_combined/
│   │   │   └── ...
│   │   ├── rf_classify/
│   │   │   └── ...
│   │   ├── process_qc/
│   │   │   └── ...
│   │   ├── combine_cta/
│   │   │   └── ...
│   │   ├── load_cta/
│   │   │   └── ...
│   │   ├── combine_clc/
│   │   │   └── ...
│   │   ├── load_clc/
│   │   │   └── ...
│   │   ├── get_meta/
│   │   │   └── ...
│   │   ├── combine_qc/
│   │   │   └── ...
│   │   └── publish_multiqc/
│   │       └── ...
│   └── nf-core/
│       └── multiqc/             # Use official nf-core MULTIQC module
│
├── bin/                         # Keep existing Python scripts
│   ├── setup.py
│   ├── get_census_adata.py
│   ├── process_query.py
│   ├── process_query_samples.py
│   ├── scvi_classify.py
│   ├── process_QC.py
│   ├── get_gemma_meta.py
│   └── utils.py
│
├── assets/
│   ├── multiqc_config.yaml      # Move from meta/
│   ├── email_template.html      # NEW: Email notification template
│   ├── sendmail_template.txt    # NEW: Plain text email
│   └── schema_input.json        # NEW: Input validation schema
│
├── docs/
│   ├── usage.md
│   ├── output.md
│   └── images/
│
├── lib/                         # NEW: Groovy helper functions
│   ├── WorkflowMain.groovy
│   └── WorkflowScannotate.groovy
│
├── .github/
│   └── workflows/
│       ├── ci.yml               # NEW: CI testing
│       ├── linting.yml          # NEW: nf-core linting
│       ├── branch.yml           # NEW: Branch protection
│       └── release.yml          # NEW: Release automation
│
├── tests/                       # Reorganized for nf-test
│   ├── nextflow.config
│   └── ...
│
├── CHANGELOG.md                 # NEW
├── CITATIONS.md                 # NEW
├── CODE_OF_CONDUCT.md           # NEW
├── LICENSE                      # MIT license
└── README.md                    # Updated
```

---

## 2. File-by-File Changes

### 2.1 Entry Point: `main.nf`

**Current:** `sc-annotate.nf` contains both processes and workflow logic (514 lines)

**Target:** Minimal `main.nf` that:
- Validates parameters using nf-validation plugin
- Includes the main workflow from `workflows/scannotate.nf`
- Handles completion summary

```nextflow
#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { SCANNOTATE } from './workflows/scannotate'
include { validateParameters; paramsHelp; paramsSummaryLog } from 'plugin/nf-validation'

// Print help message
if (params.help) {
    log.info paramsHelp("nextflow run nf-core/scannotate --input samplesheet.csv --outdir results")
    exit 0
}

// Validate input parameters
validateParameters()

workflow {
    SCANNOTATE()
}

workflow.onComplete {
    // Completion handling
}
```

### 2.2 Configuration Split

**Current:** Single `nextflow.config` (117 lines) with mixed concerns

**Target:** Split into modular configs:

| File | Purpose |
|------|---------|
| `nextflow.config` | Main config, includes others, defines profiles |
| `conf/base.config` | Default CPU/memory/time for processes |
| `conf/modules.config` | Module-specific publishDir, ext.args, etc. |
| `conf/test.config` | Minimal test data parameters |
| `conf/test_full.config` | Full-size test parameters |

### 2.3 Module Reorganization

Each process in `sc-annotate.nf` becomes a standalone module in `modules/local/TOOLNAME/`:

| Current Process | New Module Location |
|-----------------|---------------------|
| `runSetup` | `modules/local/setup_scvi/main.nf` |
| `getCensusAdata` | `modules/local/get_census_adata/main.nf` |
| `processQuery` | `modules/local/process_query_combined/main.nf` |
| `rfClassify` | `modules/local/rf_classify/main.nf` |
| `processQC` | `modules/local/process_qc/main.nf` |
| `combineCTA` | `modules/local/combine_cta/main.nf` |
| `loadCTA` | `modules/local/load_cta/main.nf` |
| `combineCLC` | `modules/local/combine_clc/main.nf` |
| `loadCLC` | `modules/local/load_clc/main.nf` |
| `getMeta` | `modules/local/get_meta/main.nf` |
| `combineQC` | `modules/local/combine_qc/main.nf` |
| `runMultiQC` | Use `modules/nf-core/multiqc/` |
| `publishMultiQC` | `modules/local/publish_multiqc/main.nf` |
| `save_params_to_file` | Remove (use nf-core standard approach) |

**Module file structure example:**
```
modules/local/rf_classify/
├── main.nf          # Process definition
├── meta.yml         # Documentation (inputs, outputs, description)
└── environment.yml  # Conda environment specification
```

### 2.4 Subworkflow Organization

Group related modules into logical subworkflows:

| Subworkflow | Modules Included |
|-------------|------------------|
| `INPUT_CHECK` | Validate samplesheet/input |
| `DOWNLOAD_STUDIES` | Existing download logic |
| `PREPARE_REFERENCE` | `setup_scvi`, `get_census_adata` |
| `PROCESS_QUERIES` | `process_query_sample`, `process_query_combined` |
| `CLASSIFY_CELLTYPES` | `rf_classify`, `combine_cta` |
| `QC_REPORTING` | `process_qc`, `combine_qc`, `multiqc` |
| `GEMMA_UPLOAD` | `load_cta`, `load_clc`, `publish_multiqc` |

### 2.5 Assets Directory

Move from `meta/` to `assets/`:

| Current Location | New Location |
|------------------|--------------|
| `meta/multiqc_config.yaml` | `assets/multiqc_config.yaml` |
| `meta/cell_type_markers.tsv` | `assets/cell_type_markers.tsv` |
| `meta/rename_cells_*.tsv` | `assets/rename_cells_*.tsv` |
| `meta/gemma_genes.tsv` | `assets/gemma_genes.tsv` |

---

## 3. New Required Files

### 3.1 `nextflow_schema.json`
- Auto-generated using `nf-core schema build`
- Defines all parameters with types, descriptions, defaults
- Powers `--help`, validation, and nf-core launch

### 3.2 `.nf-core.yml`
```yaml
repository_type: pipeline
nf_core_version: "2.14"
org_path: nf-core
lint:
  files_unchanged:
    - .github/workflows/branch.yml  # Example: skip linting certain files
```

### 3.3 `modules.json`
- Tracks installed nf-core modules and their versions
- Auto-managed by `nf-core modules install/update`

### 3.4 `CHANGELOG.md`
- Document changes per version
- Follow Keep a Changelog format

### 3.5 `CITATIONS.md`
- List all tools/references to cite
- Include scVI, CellxGene, scanpy, etc.

### 3.6 `CODE_OF_CONDUCT.md`
- Standard nf-core code of conduct

---

## 4. Container/Environment Changes

### Current Issues
- Hardcoded conda paths: `/home/rschwartz/anaconda3/envs/scanpyenv`
- Single container for all processes

### Required Changes
1. Create per-module `environment.yml` files
2. Build versioned containers for each module (or use existing biocontainers)
3. Use Seqera Containers or BioContainers where possible
4. Remove hardcoded paths

**Example `environment.yml`:**
```yaml
name: scannotate_rf_classify
channels:
  - conda-forge
  - bioconda
dependencies:
  - python=3.10
  - scanpy=1.9.6
  - scikit-learn=1.3.2
  - anndata=0.10.3
```

---

## 5. Testing Requirements

### 5.1 Test Data
- Create minimal test dataset (small h5ad files)
- Host on nf-core test-datasets or similar
- Define in `conf/test.config`

### 5.2 CI/CD Workflows
```yaml
# .github/workflows/ci.yml
- Run pipeline with test profile
- Test both conda and singularity
- Lint with nf-core tools
```

### 5.3 nf-test Integration
- Write nf-test cases for each module
- Test workflow end-to-end

---

## 6. Parameter Standardization

### Rename Parameters (nf-core conventions)
| Current | Proposed |
|---------|----------|
| `study_names` | `input` (use samplesheet) |
| `study_paths` | Remove (use samplesheet) |
| `organism` | Keep |
| `census_version` | Keep |
| `subsample_ref` | `subsample_reference` |
| `ref_collections` | `reference_collections` |
| `nmads` | `outlier_nmads` |
| `use_staging` | `gemma_staging` |

### Input Samplesheet
Replace `study_names`/`study_paths` with standardized samplesheet:
```csv
sample,study_name,study_path
sample1,GSE154208,/path/to/data
sample2,GSE123456,
```

---

## 7. Implementation Order

### Phase 1: Scaffold (1-2 days)
1. Run `nf-core create` to generate template
2. Copy existing `bin/` scripts
3. Move `meta/` contents to `assets/`
4. Create basic `conf/` structure

### Phase 2: Modules (3-5 days)
1. Extract each process into `modules/local/TOOL/main.nf`
2. Add `meta.yml` documentation for each
3. Create `environment.yml` for each module
4. Install nf-core MULTIQC module

### Phase 3: Subworkflows (2-3 days)
1. Create subworkflows grouping related modules
2. Refactor main workflow to use subworkflows
3. Implement input validation

### Phase 4: Configuration (1-2 days)
1. Split `nextflow.config` into modular configs
2. Generate `nextflow_schema.json`
3. Set up profiles (docker, singularity, conda, test)

### Phase 5: Testing & CI (2-3 days)
1. Create minimal test data
2. Write `conf/test.config`
3. Set up GitHub Actions CI
4. Add nf-test cases

### Phase 6: Documentation (1-2 days)
1. Update README.md
2. Write docs/usage.md and docs/output.md
3. Add CHANGELOG.md, CITATIONS.md
4. Generate parameter documentation

---

## 8. Key Decisions Needed

1. **Pipeline name**: Suggest `nf-core/scannotate` or `nf-core/sccellannotate`
2. **Samplesheet format**: Define required columns
3. **Container strategy**: Single container vs per-module?
4. **Test data hosting**: nf-core test-datasets or self-hosted?
5. **Gemma integration**: Keep as optional module or separate?

---

## 9. Resources

- [nf-core Pipeline Structure](https://nf-co.re/docs/contributing/pipelines/pipeline_file_structure)
- [nf-core Module Specifications](https://nf-co.re/docs/guidelines/components/modules)
- [Adding a Pipeline Tutorial](https://nf-co.re/docs/tutorials/adding_a_pipeline/overview)
- [nf-core Tools Documentation](https://nf-co.re/docs/nf-core-tools/pipelines/create)

---

## Summary

The reorganization requires:
- **Structural changes**: New directory layout following nf-core template
- **Process extraction**: Move 13+ processes to individual modules
- **Subworkflow creation**: Group modules logically
- **Configuration split**: Modular config files
- **Environment definitions**: Per-module conda/container specs
- **New files**: Schema, changelog, citations, CI workflows
- **Testing infrastructure**: Test data, nf-test, GitHub Actions