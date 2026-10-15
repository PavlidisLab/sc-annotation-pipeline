# Single-cell Re-annotation Pipeline

## Overview

This Nextflow pipeline annotates cell types in single-cell RNA-seq data from the [Gemma database](https://gemma.msl.ubc.ca) or from local files. Cells are labelled against CELLxGENE Census references, using a random forest or kNN classifier on scVI embeddings or Seurat label transfer. Assignments are made at subclass, class and family level. It also computes a QC mask, writes a MultiQC report per study and can upload the cell type assignments (CTAs), mask and report to Gemma. 

## Quick start

```bash
nextflow run main.nf -profile conda -params-file params.mm.json --input assets/samplesheet.csv
```

Use `params.hs.json` for human. Command-line arguments override the params file. The legacy `--study_names` and `--study_paths` (space-separated list or text file) replace `--input`.

The samplesheet is a CSV with the columns `sample` (required, no spaces), `study_name` (Gemma study to download) and `study_path` (MEX directory or one `.h5ad`). Each row needs `study_name` or `study_path`. A local `.h5ad` needs raw counts in `X` and Ensembl IDs in `var_names`.

## Key parameters

| Parameter | Default | Meaning |
|---|---|---|
| `classifier` | `rf` | `rf`, `knn` or `seurat` (`params.mm.json` sets `knn`) |
| `organism` | `mus_musculus` | or `homo_sapiens` |
| `census_version` | `2024-07-01` | CELLxGENE Census release |
| `subsample_ref` | `100` | Reference cells per cell type |
| `cutoff` | `0` | Minimum probability or score to assign a label |
| `ref_collections` | two mouse collections | Census collections |
| `use_gemma` | `true` | `false` skips all Gemma calls in QC and upload |
| `upload_gemma` | `true` | `false` skips only the upload |
| `upload_cta`, `upload_clc`, `upload_mask`, `upload_multiqc` | `true` | Select what to upload |
| `use_staging` | `true` | Gemma staging server |
| `cta_protocol` | `null` | Override protocol name |
| `ref_tag` | `null` | Override reference tag |

Others (`nmads`, `ref_keys`, asset files) are in `nextflow.config`.

## Classifiers

`rf` and `knn` embed query cells with the scVI model, then classify against the subsampled reference embeddings. `seurat` runs Seurat label transfer on the reference and query objects. Its R steps have no container yet, so `seurat` needs the `r4.3` conda environment (`-profile conda`). 

## Outputs

Files go to `--outdir` (default `results/<organism>_subsample_ref_<n>_<timestamp>`). Layout for `rf` and `knn` (see `conf/modules.config`):

| Path | Content |
|---|---|
| `scvi_model/`, `reference/`, `<study>/processed/` | Model, reference, query embeddings |
| `<study>/predicted_celltypes/` | CTA files |
| `<study>/masks/`, `qc/`, `combined_qc/` | QC mask, metrics |
| `multiqc/<study>/` | MultiQC report |
| `pipeline_info/` | Nextflow reports |

The `seurat` run in `results/seurat_mouse_staging_GSE124952/` instead has `ref/`, `query/`, `predict/` (scores) and `classify/` (CTA files).

## Gemma upload

Uploads use `gemma-cli`, or `gemma-cli-staging` when `use_staging` is true (the default), so set `--use_staging false` for production. Credentials come from the `GEMMA_USERNAME` and `GEMMA_PASSWORD` environment variables, never from arguments.

The CTA protocol is `sc-pipeline-<version>-<classifier>-census<census>-ref<tag>-n<subsample_ref>-cutoff<cutoff>`, e.g. `sc-pipeline-2.0.0-seurat-census20240701-ref1710b1-n100-cutoff0` (version from `params.version`). The tag is `ref_tag` or the first six hex characters of the SHA-256 of the sorted `ref_collections` (`assets/ref_tags.tsv`). The CTA name is protocol plus level. The protocol must exist in Gemma before upload:

```bash
gemma-cli addProtocol --name "sc-pipeline-2.0.0-seurat-census20240701-ref1710b1-n100-cutoff0"
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

See `markdowns/CHANGELOG.md` for changes and `markdowns/CITATIONS.md` for citations.
