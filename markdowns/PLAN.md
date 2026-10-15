# Plan for sc-annotation-pipeline

## Done

The nf-core style reorganization that this file originally planned is complete: minimal `main.nf`, `workflows/scannotate.nf`, six local subworkflows, one module per process in `modules/local/`, `conf/base.config` and `conf/modules.config`, `nextflow_schema.json`, samplesheet input (`--input`), `.nf-core.yml`, `markdowns/CHANGELOG.md`, `markdowns/CITATIONS.md` and `markdowns/CODE_OF_CONDUCT.md`, `meta.yml` for most modules, and the human/mouse test profile matrix (`conf/test_*.config`, run with `scripts/run_test_profiles.sh`). Reference files live in `assets/`. We also added the kNN classifier (`CLASSIFY_SCVI`), the `use_gemma` and `upload_gemma` switches, the Seurat label transfer classifier (`--classifier seurat`), a configurable cutoff for all three classifiers, the `subsample_ref` default of 100, and a CTA protocol built from the configuration (`sc-pipeline-<version>-<classifier>-ref<subsample_ref>-cutoff<cutoff>`, override with `--cta_protocol`).

## Removed

We deleted the original directory and file-by-file sections, the module and subworkflow mapping tables, the phase schedule, and the key decisions list, because the code now has the structure they proposed. The proposed parameter renames (`input`-only, `outlier_nmads`, `gemma_staging`) were dropped: the current names are in use in all params files.

## Open

1. Release v2.1.0 from `development` (branch `release-2.1.0`), then v2.2.0 with the Seurat classifier from `seurat-classifier`.
2. Protocol nomenclature. The current identifier omits the census version and the reference collections. A fuller scheme with the census version and a reference hash tag is planned, not implemented.
3. Per-reference building, to match the benchmark's `split_refs` (one reference per collection instead of a pooled reference).
4. Missing ontology URIs in `assets/rename_cells_mus_musculus.tsv` for some mouse types (the `deep layer non-IT` class and several hippocampal subclasses).
5. No CI and no tests. A local `Jenkinsfile` (gitignored) is dormant, and there is no `tests/` directory, nf-test cases, or GitHub Actions workflow. Run `scripts/run_test_profiles.sh` manually before a release.
6. Per-module `environment.yml` files and containers. Modules still point at hardcoded conda paths under `/home/rschwartz/anaconda3/envs/`.
