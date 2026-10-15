# Single-cell Re-annotation Pipeline

## Abstract

This Nextflow pipeline (version 2.1.0) annotates cell types in single-cell RNA-seq data from the [Gemma database](https://gemma.msl.ubc.ca) or from local files. Cells are labelled against CELLxGENE Census references with a random forest (`rf`) or kNN (`knn`) classifier on scVI embeddings. The pipeline also computes a QC mask and a MultiQC report per study, and can upload the results to Gemma.

![Workflow of the scannotate pipeline](docs/diagrams/scannotate_flow.png)

*Figure 1. Workflow of the pipeline (source: `docs/diagrams/scannotate_flow.mmd`). Blue nodes are inputs, green nodes are pipeline steps and amber nodes are outputs; arrows show data flow. The figure shows how studies and the Census reference are embedded with the scVI model, classified, scored for QC and optionally uploaded to Gemma.*

## Methods

**Inputs.** We read studies from a samplesheet, from Gemma (`study_name`) or from local paths (`study_path`). A local `.h5ad` needs raw counts in `X` and Ensembl IDs in `var_names`. Cells with zero counts are removed and logged (`bin/process_query.py`).

**Reference preparation.** We drew `subsample_ref` cells per cell type (default 500) from the Census release `census_version` (default `2024-07-01`) and the collections in `ref_collections`, and embedded studies and reference with the scVI model.

**Classification.** We trained `rf` (default) or `knn` on the reference and assigned labels at subclass, class and family level. A label was assigned when its probability reached `cutoff` (default 0). Defaults are in `nextflow.config`.

**QC.** We computed a QC mask with median absolute deviation thresholds (`nmads`) and wrote one MultiQC report per study.

**Gemma upload.** Upload is off unless `use_gemma` and `upload_gemma` are both true. Uploads use `gemma-cli-staging` when `use_staging` is true (the default). The CTA protocol is `sc-pipeline-2.1.0` and the CTA name is `sc-pipeline-2.1.0-<level>`. Each upload replaces the CTA of the same name. `upload_cta`, `upload_clc`, `upload_mask` and `upload_multiqc` select what is sent.

## Usage

**Prerequisites.** Install Nextflow 23.04.0 or later (`manifest.nextflowVersion`). The default `conda` profile builds environments itself; `singularity`, `apptainer` and `docker` profiles are also available. The default executor is Slurm (`process.executor`). Gemma credentials are needed only when `use_gemma` is true.

1. Clone the repository and choose a profile.

```bash
git clone https://github.com/PavlidisLab/sc-annotation-pipeline.git
cd sc-annotation-pipeline && git checkout release-2.1.0
```

2. To run the local test, we used the `test_mouse_local` profile. It reads one local MEX directory and sets `--use_gemma false`, so it needs no credentials.

```bash
nextflow run main.nf -profile conda,test_mouse_local
```

3. To run real data, write a samplesheet with the columns `sample` (required, no spaces), `study_name` (Gemma study to download) and `study_path` (MEX directory or one `.h5ad`). Each row needs `study_name` or `study_path`. The example row from `assets/samplesheet.csv` downloads study GSE124952:

```
sample,study_name,study_path
GSE124952,GSE124952,
```

```bash
nextflow run main.nf -profile conda --input assets/samplesheet.csv \
  --organism mus_musculus --census_version 2024-07-01 --outdir results/my_run
```

`--organism` is `mus_musculus` or `homo_sapiens`, `--census_version` is the Census release and `--outdir` is where outputs go. The schema marks none as required, and the defaults come from `nextflow.config`. Add `-params-file params.mm.json` (or `params.hs.json`) for preset values; command-line flags override it.

4. To upload to Gemma, follow this order. Create the protocol, set credentials, then run.

```bash
gemma-cli-staging addProtocol --name "sc-pipeline-2.1.0"
export GEMMA_USERNAME=<user> GEMMA_PASSWORD=<password>
nextflow run main.nf -profile singularity --input assets/samplesheet.csv --use_gemma true --upload_gemma true
```

Pass `--use_staging false` for production and use `gemma-cli` for the protocol. Credentials come from the environment only.

**Troubleshooting.**

- `required environment variable(s) not set: GEMMA_USERNAME, GEMMA_PASSWORD` means the variables are unset when `use_gemma` is true. Export them or set `--use_gemma false`.
- With `singularity` or `apptainer`, `envWhitelist` passes the credentials into the container by name, so export them before the run.
- MultiQC is absent when `process_samples` is true; set it to false to get reports.
- A sample with zero counts is removed and logged to stderr; the run stops only if every cell is empty.

## Outputs

Files go to `--outdir` (default `results/<organism>_subsample_ref_<n>_<timestamp>`), laid out by `conf/modules.config`. A real run (`results/mus_musculus_subsample_ref_500_2026-10-15_13-21-31`) held these folders:

- `download/`, `get/`: downloaded study data and `<study>_sample_meta.tsv`.
- `scvi_model/`, `reference/`: scVI model and Census reference.
- `<study>/processed/`, `predicted_celltypes/`, `qc/`: query embeddings, CTA files, QC metrics and mask.
- `multiqc/<study>/`: MultiQC report.
- `pipeline_info/`: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt`, `pipeline_dag.dot`.

## Limitations

There is no CI or automated test suite. `scripts/run_test_profiles.sh` runs the eight test profiles with Singularity and `--upload_gemma false`. `GET_META` runs on the launch host, and credential inheritance by Slurm tasks was not tested. MultiQC is produced only when `process_samples` is false.

## References

See `markdowns/CHANGELOG.md` for changes and `markdowns/CITATIONS.md` for citations.
