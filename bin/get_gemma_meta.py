#!/bin/python
import warnings
warnings.filterwarnings("ignore")
import pandas as pd
import argparse
import gemmapy

def argument_parser():
    parser = argparse.ArgumentParser(description="Preprocess data from GEMMA")
    parser.add_argument("--study_name", type=str, help="Name of the study", default="GSE152715.1")
    parser.add_argument('--gemma_username', type=str, default=None)
    parser.add_argument('--gemma_password', type=str, default=None)
    return parser.parse_args()

def main():
    args = argument_parser()

    gemma_username = args.gemma_username
    gemma_password = args.gemma_password
    client = gemmapy.GemmaPy(auth=[gemma_username,gemma_password], path='staging')
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

    rows = []
    for s in samples_raw.data:
        row = {
            "sample_id": s.id,             # BioAssay ID (join key)
            "biomaterial_id": s.sample.id,  # BioMaterial ID (reference only)
            "sample_name": s.name,
            "organism": s.array_design.taxon.scientific_name.lower().replace(" ", "_"),
        }
        for c in s.sample.characteristics:
            if c.category is not None:
                row[c.category] = c.value
        rows.append(row)

    sample_meta_df = pd.DataFrame(rows)
    sample_meta_df.to_csv(f"{study_name}_sample_meta.tsv", index=False, sep='\t')
    

    
if __name__ == "__main__":
    main()
