# Copilot Instructions for cell_annotation_cortex.nf

## Project Overview
This repository contains a Nextflow pipeline for automated cell type annotation of single-cell transcriptomics data, integrating reference datasets and models from the CellxGene census and Gemma database. The pipeline leverages scVI embeddings and a random forest classifier for label transfer, with extensive QC and reporting.

## Architecture & Data Flow
- **Pipeline Orchestration:** Main workflow in `sc-annotate.nf` (see also `nextflow.config`).
- **Python Modules:** All core data processing, model inference, and QC are handled by scripts in `bin/` (notably `scvi_classify.py`, `process_QC.py`, `get_census_adata.py`, etc.).
- **Reference Data:** Downloaded and processed via CellxGene census APIs; reference collections are specified in `params.*.json`.
- **Outputs:** Results are written to `results/` with per-run subfolders, including predicted cell types, QC reports, and parameter logs.

## Key Workflows
- **Run the pipeline:**
  ```sh
  nextflow run sc-annotate.nf -profile conda -params-file params.mm.json --study_names "GSE12345 GSE67890"
  ```
  Or use `--study_paths` for local MEX files. See `README.md` for more details.
- **Resume after error:**
  ```sh
  nextflow run sc-annotate.nf -profile conda -resume -params-file params.mm.json
  ```
- **Parameterization:**
  - Use `params.mm.json` (mouse) or `params.hs.json` (human) as templates. Do not edit these directly; copy and modify as needed.
  - Command-line arguments override JSON and `nextflow.config` defaults.

## Project Conventions
- **Environment:** Uses Conda environments (see `nextflow.config`). Some paths are hard-coded; update as needed for new deployments.
- **Reference Collections:** Always specify via JSON due to Nextflow CLI parsing limitations with spaces.
- **QC & Reporting:** MultiQC reports and outlier detection are handled in `bin/process_QC.py` using MAD-based thresholds.
- **Cell Type Mapping:** Taxonomy mapping files are in `meta/` and referenced via parameters.
- **Results:** Each run creates a timestamped subfolder in `results/`.

## Integration Points
- **External APIs:** CellxGene census (for reference data/models), Gemma CLI (for data upload).
- **Containerization:** Some processes specify Docker/Singularity containers in `nextflow.config`.

## Examples & References
- See `README.md` for parameter explanations, input/output structure, and example commands.
- Example output: `results/<run>/GSE*/predicted_celltypes/`.
- Example MultiQC: `images/multiqc/`.

## Tips for AI Agents
- Always check `params.*.json` and `nextflow.config` for current defaults and required parameters.
- When adding new processes, follow the pattern in `sc-annotate.nf` and use Python scripts in `bin/`.
- For new reference data, update the relevant JSON and mapping files in `meta/`.
- For debugging, inspect per-run logs and parameter files in the output directory.

---
For further details, see the full [README.md](../README.md).
