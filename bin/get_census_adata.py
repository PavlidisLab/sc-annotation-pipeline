#!/user/bin/python3

import warnings
warnings.filterwarnings("ignore")
from pathlib import Path
import random
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
from sklearn.ensemble import RandomForestClassifier
import adata_functions
from adata_functions import *
from pathlib import Path
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, classification_report, confusion_matrix, roc_curve, auc
from sklearn.preprocessing import label_binarize
import matplotlib.pyplot as plt
import seaborn as sns
import json
import argparse
import os
import json

# Function to parse command line arguments
def parse_arguments():
    parser = argparse.ArgumentParser(description="Download model file based on organism, census version, and tree file.")
    parser.add_argument('--organism', type=str, default='mus_musculus', help='Organism name (e.g., homo_sapiens)')
    parser.add_argument('--census_version', type=str, default='2024-07-01', help='Census version (e.g., 2024-07-01)')
    parser.add_argument('--ref_collections', type=str, nargs = '+', default = [
        "A taxonomy of transcriptomic cell types across the isocortex and hippocampal formation",
        "An integrated transcriptomic and epigenomic atlas of mouse primary motor cortex cell types",
        "Adult mouse cortical cell taxonomy revealed by single cell transcriptomics",
        "Tabula Muris Senis",
        "Single-cell transcriptomics characterization of oligodendrocytes and microglia in white matter aging"
    ]) 
    parser.add_argument('--seed', type=int, default=42)
    parser.add_argument('--organ', type=str, default="brain")
    parser.add_argument('--assay', type=str, nargs = "+", help="Assays to subset from referenc (unnecessary)", default=None)
    parser.add_argument('--tissue', type=str, nargs="+", default = None, help = "tissues to pull from (different from organ, this can select for more specific brain regions)")
    parser.add_argument('--subsample', type=int, help="Number of cells per cell type to subsample from reference", default=50)
    parser.add_argument('--rename_file', type=str, default="/space/grp/rschwartz/rschwartz/cell_annotation_cortex.nf/meta/author_cell_annotations/rename_cells_mmus_author.tsv")
    parser.add_argument('--ref_name', type=str, default="whole_cortex", help="Prefix of temporary reference file created")
    parser.add_argument('--original_celltype_columns', type=str, default=None)
    parser.add_argument('--author_annotations_path', type=str, default=None)
    
    if __name__ == "__main__":
        known_args, _ = parser.parse_known_args()
        return known_args

def main():
   # Parse command line arguments
   args = parse_arguments()


   # Set organism and census_version from arguments
   organism = args.organism
   census_version = args.census_version
   ref_collections=args.ref_collections
   SEED = args.seed
   assay = args.assay
   subsample = args.subsample
   tissue = args.tissue
   organ=args.organ
   rename_file=args.rename_file
   original_celltype_columns = args.original_celltype_columns
   author_annotations_path = args.author_annotations_path
   
   #f organism == "mus_musculus":
   if original_celltype_columns is not None:
      original_celltypes = get_original_celltypes(columns_file=original_celltype_columns,
                                          author_annotations_path=author_annotations_path) 
   else:
      original_celltypes = None
  
   ref=get_census(organism=organism, 
                     subsample=subsample, census_version=census_version, organ=organ,
                        ref_collections=ref_collections, assay=assay, tissue=tissue, rename_file=rename_file, seed=SEED, original_celltypes=original_celltypes)

   print("finished fetching anndata")
   outdir="refs"
   os.makedirs(outdir, exist_ok=True) 

    # handle tissue and assay params, can be lists
   if tissue and assay:
      ref_name = "_".join([f"{t}-{a}" for t in tissue for a in assay])
   else:
      ref_name = args.ref_name

   if len(ref.obs.index) == 0:
      raise ValueError(f"Reference {ref_name} has no cells, check README for proper ref collections")
   if len(ref.var.index) == 0:
      raise ValueError(f"Reference {ref_name} has no genes, check README for proper ref collections")
   else:
      new_ref_name = ref_name.replace(" ", "_").replace("\\/", "_") \
      .replace("(","").replace(")","") \
      .replace("\\", "") \
      .replace("'", "") \
      .replace(":", "")
      ref.write(os.path.join(outdir,f"{new_ref_name}.h5ad"))

if __name__ == "__main__":
    main()