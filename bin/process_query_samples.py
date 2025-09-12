#!/user/bin/python3

from pathlib import Path
import os
import sys
import scanpy as sc
import numpy as np
import pandas as pd
import anndata as ad
import cellxgene_census
from scipy.sparse import csr_matrix
import warnings
import cellxgene_census
import cellxgene_census.experimental
import scvi
from utils import *
from pathlib import Path
import argparse
import os
import json
import scipy
import gzip
    
def parse_arguments():
  parser = argparse.ArgumentParser(description="Download model file based on organism, census version, and tree file.")
  parser.add_argument('--model_path', type=str, default="/space/grp/rschwartz/rschwartz/cell_annotation_cortex.nf/work/31/eda2bebc822637d00e3b8963b1d11f/scvi-homo_sapiens-2024-07-01", help='Path to the scvi model file')
  parser.add_argument('--query_name', type=str, default="1011640_GSM7050280", help='Name of the study')
  parser.add_argument('--query_path', type=str, default="/space/grp/rschwartz/rschwartz/cell_annotation_cortex.nf/work/5a/8b93c45f0cf99abbc07813af44d247/GSE225554/1011640_GSM7050280", help='Path to the study file')
  parser.add_argument('--seed', type=int, default=42)
   
  if __name__ == "__main__":
    known_args, _ = parser.parse_known_args()
    return known_args


def has_expression_data(adata):
    return hasattr(adata, "X") and getattr(adata.X, "nnz", None) not in (0, None)

def check_size(adata):
  if adata.shape[0] < 50:
    return False
  else:
    return True
  

def main():

  args = parse_arguments()
  query_path = args.query_path
  query_name = args.query_name
  sample_id = query_name.split("_")[0]
  model_path = args.model_path
  SEED = args.seed
  random.seed(SEED)         # For `random`
  np.random.seed(SEED)      # For `numpy`
  scvi.settings.seed = SEED # For `scvi`
  
  
  try:
    # Attempt to read the 10x mtx data
    adata = sc.read_10x_mtx(query_path)
    adata.obs["sample_id"] = sample_id  # Add sample_id to obs
    adata.obs_names_make_unique()
  except Exception as e:
    print(f"Error processing {sample_id} automatically: {e}. Trying manual read.")
    
    # If an error occurs, try reading the files manually
    try:
        # Read the matrix, genes, and barcodes files manually
        matrix_path = os.path.join(query_path, "matrix.mtx.gz")
        genes_path = os.path.join(query_path, "features.tsv.gz")
        barcodes_path = os.path.join(query_path, "barcodes.tsv.gz")
        
        # Load the matrix in CSR format
        with gzip.open(matrix_path, 'rb') as f:
            matrix = scipy.io.mmread(f).tocsr()
        
        # Read the gene and barcode files
        with gzip.open(genes_path, 'rt') as f:
            genes = [line.strip().split("\t") for line in f]
        with gzip.open(barcodes_path, 'rt') as f:
            barcodes = [line.strip() for line in f]
        
        # Create AnnData object
        adata = sc.AnnData(X=matrix.T)  # Transpose to match expected shape (cells x genes)
        adata.var_names = [gene[1] for gene in genes]
        adata.var_names_make_unique()# gene ids as the variable names
        adata.obs_names = barcodes  # cell barcodes as the observation names
        adata.obs["sample_id"] = sample_id  # Add sample_id to obs
        adata.obs_names_make_unique()  # make sure the observation names are unique
        # Store the AnnData object in the dictionary
     #   all_sample_ids[new_sample_id] = adata
        print(f"Successfully created AnnData for {sample_id} from individual files.")
    
    except Exception as manual_e:
      adata = sc.AnnData(X=csr_matrix((0, 0)))
      adata.obs["sample_id"] = sample_id  # Add sample_id to obs
      adata.write_h5ad(f"{sample_id}_empty.h5ad")
      #print(f"Error processing {sample_id} manually: {manual_e}")
      raise Warning(f"Failed to process {sample_id} using both automatic and manual methods: {manual_e}")

  if has_expression_data(adata) is False:
    adata = sc.AnnData(X=csr_matrix((0, 0)))
    # write a fake h5ad to trick nextflow
    adata.write_h5ad(f"{sample_id}_empty.h5ad")
    raise Warning(f"Sample {sample_id} has no expression data. Skipping.")
  if check_size(adata) is False:
    os.makedirs("small_samples", exist_ok=True)
    adata.obs["cell_id"] = adata.obs.index
    adata.write_h5ad(os.path.join("small_samples",f"{query_name}.h5ad"))
  else:
    adata.obs["cell_id"] = adata.obs.index
    adata.write_h5ad(f"{query_name}_raw.h5ad")
    adata = process_query(adata, model_path, batch_key="sample_id", seed=SEED)
    adata.write_h5ad(f"{query_name}.h5ad")
      
        
if __name__ == "__main__":
    main()
    