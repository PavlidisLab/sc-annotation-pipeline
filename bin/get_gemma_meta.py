#!/bin/python
import warnings
warnings.filterwarnings("ignore")
import pandas as pd
import anndata as ad
from scipy.sparse import csr_matrix
from pathlib import Path
import argparse
import os
import json
import sys
import numpy as np
import gemmapy

def argument_parser():
    parser = argparse.ArgumentParser(description="Preprocess data from GEMMA")
    parser.add_argument("--study_name", type=str, help="Name of the study", default="GSE152715.1")
    parser.add_argument('--gemma_username', type=str, default="raschwar")
    parser.add_argument('--gemma_password', type=str, default="7nddtt")
    return parser.parse_args()

def main():
    args = argument_parser()

    gemma_username = args.gemma_username
    gemma_password = args.gemma_password
    client = gemmapy.GemmaPy(auth=[gemma_username,gemma_password], path='dev')
    study_name = args.study_name
    samples_raw = client.raw.get_dataset_samples(study_name)
  
    samples = client.get_dataset_samples(study_name, use_processed_quantitation_type=False)
    sample_names = [x for x in samples["sample_name"]]
    sample_ids = [x.id for x in samples_raw.data]
    organisms = [x.array_design.taxon.scientific_name.lower().replace(" ", "_") for x in samples_raw.data]
    sample_meta = [df for df in samples["sample_characteristics"]]
    # add sample id to each df in sample meta
    sample_meta_updated1 = [df.assign(sample_id=sample_ids[i]) for i, df in enumerate(sample_meta)]
    # add organism to each df in sample meta
    sample_meta_updated2 = [df.assign(organism=organisms[i]) for i, df in enumerate(sample_meta_updated1)]
    # add names back
    sample_meta_updated2 = [df.assign(sample_name=sample_names[i]) for i, df in enumerate(sample_meta_updated2)]
    
    # combine all dfs
    sample_meta_combined = pd.concat(sample_meta_updated2)
    sample_meta_combined.drop_duplicates(subset=["sample_id", "category"], inplace=True)
    sample_meta_df = sample_meta_combined.pivot(index=["sample_id","sample_name","organism"], columns="category",values="value").reset_index()
    sample_meta_df.to_csv(f"{study_name}_sample_meta.tsv", index=False, sep='\t')
    

    
if __name__ == "__main__":
    main()
