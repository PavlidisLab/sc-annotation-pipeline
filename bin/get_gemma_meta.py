#!/bin/python
import warnings
warnings.filterwarnings("ignore")
import pandas as pd
import argparse
import os
import sys

def argument_parser():
    parser = argparse.ArgumentParser(description="Preprocess data from GEMMA")
    parser.add_argument("--study_name", type=str, help="Name of the study", default="GSE152715.1")
    parser.add_argument('--use_staging', action='store_true', help='Use Gemma staging instead of production')
    return parser.parse_args()

def main():
    args = argument_parser()

    # Credentials come from the environment only, so they never appear on a command line
    # or in work/*/.command.sh. Fail before any network call if either is missing.
    missing = [v for v in ("GEMMA_USERNAME", "GEMMA_PASSWORD") if not os.environ.get(v)]
    if missing:
        sys.exit(f"ERROR: required environment variable(s) not set: {', '.join(missing)}")
    gemma_username = os.environ["GEMMA_USERNAME"]
    gemma_password = os.environ["GEMMA_PASSWORD"]
    # gemmapy treats path as a host alias, so production is the default (no path)
    import gemmapy  # imported after the credential check so a missing variable fails first
    client_kwargs = {'path': 'staging'} if args.use_staging else {}
    client = gemmapy.GemmaPy(auth=[gemma_username,gemma_password], **client_kwargs)
    study_name = args.study_name

    # use_processed_quantitation_type=False is the raw per-BioAssay set this
    # pipeline needs (matches MEX/h5ad filenames, CTA files); each object
    # already carries its own BioMaterial (`.sample`) with characteristics
    # nested inside, so id/name/characteristics all come from one place with
    # no risk of misalignment. (=True instead explodes each sample into one
    # BioAssay per cell type, unrelated to what's needed here.) A previous
    # version of this script instead fetched names/characteristics from a
    # *second*, separately-ordered API call and zipped the two together by
    # list position; that second call's own id column turned out to be a
    # BioMaterial ID (a different ID space than BioAssay, zero overlap), so
    # the positional zip silently attached the wrong sample's name/
    # characteristics to a given BioAssay id in several studies.
    samples_raw = client.raw.get_dataset_samples(study_name, use_processed_quantitation_type=False)

    # generic platforms (e.g. Generic_mouse_ncbiIds) carry no taxon, so fall back to the dataset's own taxon
    dataset_taxon = client.raw.get_datasets_by_ids([study_name]).data[0].taxon.scientific_name

    rows = []
    for s in samples_raw.data:
        taxon = s.array_design.taxon.scientific_name if s.array_design.taxon is not None else dataset_taxon
        row = {
            "sample_id": s.id,             # BioAssay ID (join key)
            "biomaterial_id": s.sample.id,  # BioMaterial ID (reference only)
            "sample_name": s.name,
            "organism": taxon.lower().replace(" ", "_"),
        }
        for c in s.sample.characteristics:
            if c.category is not None:
                row[c.category] = c.value
        rows.append(row)

    sample_meta_df = pd.DataFrame(rows)
    sample_meta_df.to_csv(f"{study_name}_sample_meta.tsv", index=False, sep='\t')
    

    
if __name__ == "__main__":
    main()
