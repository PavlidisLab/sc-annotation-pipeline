# Single-cell Re-annotation Pipeline

## Overview

This Nextflow pipeline annotates cell types in single-cell RNA-seq data from the [Gemma database](https://gemma.msl.ubc.ca) or from local files. Cells are labelled against CELLxGENE Census references with a random forest or kNN classifier on scVI embeddings. Assignments are made at subclass, class and family level. The pipeline also computes a QC mask, writes a MultiQC report per study and can upload the cell type assignments (CTAs), mask and report to Gemma.

## Quick start

```bash
nextflow run main.nf -profile conda -params-file params.mm.json --input assets/samplesheet.csv
```

Use `params.hs.json` for human. Command-line arguments override the params file. The legacy `--study_names` and `--study_paths` (space-separated list or text file) replace `--input`.

The samplesheet is a CSV with the columns `sample` (required, no spaces), `study_name` (Gemma study to download) and `study_path` (MEX directory or one `.h5ad`). Each row needs `study_name` or `study_path`. A local `.h5ad` needs raw counts in `X` and Ensembl IDs in `var_names`.

## Key parameters

| Parameter | Default | Meaning |
|---|---|---|
| `classifier` | `rf` | `rf` or `knn` (`params.mm.json` sets `knn`) |
| `organism` | `mus_musculus` | or `homo_sapiens` |
| `census_version` | `2024-07-01` | CELLxGENE Census release |
| `subsample_ref` | `500` | Reference cells per cell type (`params.mm.json` sets 100) |
| `cutoff` | `0` | Minimum probability to assign a label |
| `process_samples` | `false` | Process samples individually |
| `ref_collections` | two mouse collections | Census collections |
| `use_gemma` | `true` | `false` skips all Gemma calls in QC and upload |
| `upload_gemma` | `true` | `false` skips only the upload |
| `upload_cta`, `upload_clc`, `upload_mask`, `upload_multiqc` | `true` | Select what to upload |
| `use_staging` | `true` | Gemma staging server |
| `preferredCtaLevel` | `class` | Level flagged as the preferred CTA |

Others (`nmads`, `ref_keys`, `seed`, asset files) are in `nextflow.config`.

## Outputs

Files go to `--outdir` (default `results/<organism>_subsample_ref_<n>_<timestamp>`). The layout follows `conf/modules.config`:

| Path | Content |
|---|---|
| `scvi_model/`, `reference/`, `<study>/processed/` | Model, reference, query embeddings |
| `<study>/predicted_celltypes/` | CTA files |
| `<study>/masks/`, `qc/`, `combined_qc/` | QC mask, metrics |
| `multiqc/<study>/` | MultiQC report |
| `pipeline_info/` | Nextflow reports |

MultiQC is only produced when `process_samples` is false.

## Gemma upload

Uploads use `gemma-cli`, or `gemma-cli-staging` when `use_staging` is true (the default), so set `--use_staging false` for production. Credentials come from the `GEMMA_USERNAME` and `GEMMA_PASSWORD` environment variables, never from arguments. The `singularity` and `apptainer` profiles pass these variables into containers by name (`envWhitelist`), and `GET_META` runs on the launch host (`conf/modules.config`) because SLURM jobs do not inherit them.

The CTA protocol is `sc-pipeline-<version>` and the CTA name is `sc-pipeline-<version>-<level>`, with the version taken from `params.version` (2.1.0). Each upload replaces the CTA of the same name. The protocol must exist in Gemma before upload:

```bash
gemma-cli addProtocol --name "sc-pipeline-2.1.0"
```

## Profiles and testing

`conda` is on by default. `singularity`, `docker` and `apptainer` also exist. SLURM is the default executor.

| Profile | Test |
|---|---|
| `test_mouse`, `test_human` | Samplesheet input |
| `test_mouse_persample`, `test_human_persample` | `--process_samples true` |
| `test_mouse_local`, `test_human_local` | `--study_paths`, `--use_gemma false` |
| `test_mouse_upload_off`, `test_human_upload_off` | `--upload_gemma false` |

`scripts/run_test_profiles.sh` runs all eight with Singularity and `--upload_gemma false`, then reports pass or fail. There is no CI and no automated test suite.

## Versions and citation

See `markdowns/CHANGELOG.md` for changes and `markdowns/CITATIONS.md` for citations. Text removed from this README is in the git history.
