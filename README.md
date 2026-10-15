# Single-cell Re-annotation Pipeline

## Abstract

This Nextflow pipeline (version 2.1.0) annotates cell types in single-cell RNA-seq data from the [Gemma database](https://gemma.msl.ubc.ca) or from local files. Cells are labelled against CELLxGENE Census references with a random forest or kNN classifier on scVI embeddings. The pipeline also computes a QC mask and a MultiQC report per study, and can upload the results to Gemma.

![Workflow of the scannotate pipeline](docs/diagrams/scannotate_flow.png)

*Figure 1. Workflow of the pipeline (source: docs/diagrams/scannotate_flow.mmd). Blue nodes are inputs, green nodes are steps and amber nodes are outputs; arrows show data flow. Studies and the Census reference are embedded with the scVI model, classified, scored for QC and optionally uploaded to Gemma.*

## Methods

The pipeline ran six stages.

1. Input validation. Studies were read from a samplesheet as Gemma study names or local paths, and downloaded when needed.
2. Reference preparation. We downloaded the scVI model for the organism and drew reference cells from CELLxGENE Census.
3. Query processing. We embedded each query study with the scVI model.
4. Cell type classification. A random forest (default) or kNN classifier trained on the reference embeddings labelled each query cell at subclass, class and family level.
5. QC reporting. We computed outlier metrics and a QC mask and wrote one MultiQC report per study.
6. Gemma upload (optional). Cell type assignments, QC outputs and the report were uploaded to Gemma.

Cells with no counts were removed before stages 3 to 5. Reference size, Census release and cutoff are parameters (see the parameter schema for defaults).

## Usage

Install Nextflow 23.04.0 or later. The default executor is Slurm. Profiles for conda (default), singularity, apptainer and docker are available. Gemma credentials are needed only for uploads.

Nextflow runs the pipeline directly from GitHub, with no clone or checkout. The default branch (main) is kept at the latest release; adding -r with a release tag (for example -r v2.1.0) pins a version for reproducibility.

A samplesheet has the columns sample, study_name (a Gemma study) and study_path (a MEX directory or h5ad file); each row needs one of the two. This example downloads study GSE124952:

```bash
nextflow run PavlidisLab/sc-annotation-pipeline -profile conda \
  --input https://raw.githubusercontent.com/PavlidisLab/sc-annotation-pipeline/main/assets/samplesheet.csv \
  --organism mus_musculus --census_version 2024-07-01 --outdir results/my_run
```

For your own data, give --input the full path of your samplesheet. Presets for mouse and human are in params.mm.json and params.hs.json; after nextflow pull PavlidisLab/sc-annotation-pipeline, add one (flags override it):

```bash
-params-file ~/.nextflow/assets/PavlidisLab/sc-annotation-pipeline/params.mm.json
```

The test profiles use repository files, stored in ~/.nextflow/assets/PavlidisLab/sc-annotation-pipeline after the first run. The mouse local test profile reads a fixed path on the Pavlidis lab file system, so it runs only there and needs no credentials:

```bash
nextflow run PavlidisLab/sc-annotation-pipeline -profile conda,test_mouse_local
```

To upload, create the protocol, set credentials in the environment, then run. For production add --use_staging false and use gemma-cli.

```bash
gemma-cli-staging addProtocol --name "sc-pipeline-2.1.0"
export GEMMA_USERNAME=<user> GEMMA_PASSWORD=<password>
nextflow run PavlidisLab/sc-annotation-pipeline -profile singularity \
  --input /full/path/to/samplesheet.csv --use_gemma true --upload_gemma true
```

If a run stops because the Gemma variables are unset, export them or set --use_gemma false. MultiQC reports are produced only when per-sample processing is off.

## Outputs

Files go to the output directory. A run held downloaded data, the scVI model and Census reference, and per study the embeddings, CTA files, QC metrics and mask, and the MultiQC report. A pipeline_info folder holds the Nextflow execution report, timeline, trace and DAG.

## Limitations

There is no automated test suite. The script scripts/run_test_profiles.sh runs the eight test profiles with uploads disabled. The metadata step runs on the launch host, and credential inheritance by Slurm tasks was not tested.

## References

See markdowns/CHANGELOG.md for changes and markdowns/CITATIONS.md for citations.
