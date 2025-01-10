#!/user/bin/python3

from pathlib import Path
import os
import sys
import scanpy as sc
import numpy as np
import pandas as pd
import anndata as ad
import cellxgene_census
import scvi
from scipy.sparse import csr_matrix
import warnings
import cellxgene_census
import cellxgene_census.experimental
import scvi
import adata_functions
from adata_functions import *
from pathlib import Path
import argparse
import os
import json

# Function to parse command line arguments
def parse_arguments():
  parser = argparse.ArgumentParser(description="Download model file based on organism, census version, and tree file.")
  parser.add_argument('--model_path', type=str, default="/space/grp/rschwartz/rschwartz/biof501_proj/scvi-human-2024-07-01", help='Path to the scvi model file')
  parser.add_argument('--study_name', type=str, default="GSE198014", help='Name of the study')
  parser.add_argument('--study_path', type=str, default="/space/scratch/gemma-single-cell-data-ensembl-id/GSE198014", help='Path to the study file')
  parser.add_argument('--seed', type=int, default=42)
   
  if __name__ == "__main__":
    known_args, _ = parser.parse_known_args()
    return known_args

def main():

  # Parse command line arguments
  args = parse_arguments()
  SEED = args.seed
  random.seed(SEED)         # For `random`
  np.random.seed(SEED)      # For `numpy`
  scvi.settings.seed = SEED # For `scvi`
  # Set organism and census_version from arguments
  model_path = args.model_path
  study_name = args.study_name
  study_path = args.study_path
  scvi.settings.seed = args.seed # For `scvi`

  sample_ids = os.listdir(study_path)
  
  all_sample_ids = {}
  
  for sample_id in sample_ids:
    query_path = os.path.join(study_path, sample_id)
    all_sample_ids[sample_id] = sc.read_10x_mtx(query_path)
    
  #for sample_id_name, sample_id in all_sample_ids.items():
    #  sc.pp.calculate_qc_metrics(sample_id. qc_vars=["mt", "ribo", "hb"], )
    #  sc.pp.filter_genes()
    #  sc.pp.filter_cells(sample_id, min_genes=200)
    
  
  combined_adata = sc.concat(all_sample_ids, label="sample_id", join="inner")
  combined_adata.obs["cell_id"] = combined_adata.obs.index
  combined_adata.obs_names = combined_adata.obs["cell_id"].astype(str) + "_" + combined_adata.obs["sample_id"].astype(str)

  combined_adata = process_query(combined_adata, model_path, batch_key="sample_id")
  combined_adata.write_h5ad(f"{study_name}.h5ad")
  
  #sc.pp.neighbors(combined_adata, use_rep="scvi")
  #sc.tl.umap(combined_adata)
  #sc.pl.umap(combined_adata, color=["sample_id"])
  
  
  
if __name__ == "__main__":
    main()