#!/user/bin/python3

from pathlib import Path
import os
import sys
import scanpy as sc
import numpy as np
import pandas as pd
import anndata as ad
import re
from scipy.sparse import csr_matrix
import warnings
from matplotlib.patches import Patch
from collections import defaultdict
import matplotlib.pyplot as plt
import seaborn as sns
import json
import argparse
import os
import json
from types import SimpleNamespace
from scipy.stats import median_abs_deviation
import subprocess
import seaborn as sns
import scipy.spatial.distance as dist
import scipy.cluster.hierarchy as hierarchy
import adata_functions
from adata_functions import *

# Function to parse command line arguments
def parse_arguments():
    parser = argparse.ArgumentParser(description="Classify cells given 1 ref and 1 query")
    parser.add_argument('--organism', type=str, default='mus_musculus', help='Organism name (e.g., homo_sapiens)')
    parser.add_argument('--query_path', type=str, default="/space/grp/rschwartz/rschwartz/cell_annotation_cortex.nf/mmus/2a/69d63b6c1090d6f10a71cc5662301c/GSE152715.1.h5ad")
    parser.add_argument('--assigned_celltypes_path', type=str, default="/space/grp/rschwartz/rschwartz/cell_annotation_cortex.nf/mmus/2a/69d63b6c1090d6f10a71cc5662301c/GSE152715.1_predicted_celltype.tsv")
    parser.add_argument('--markers_file', type=str, default="/space/grp/rschwartz/rschwartz/cell_annotation_cortex.nf/meta/cell_type_markers.tsv")
    parser.add_argument('--gene_mapping', type=str, default="/space/grp/rschwartz/rschwartz/cell_annotation_cortex.nf/meta/gemma_genes.tsv")
    parser.add_argument('--nmads',type=int, default=5)
    parser.add_argument('--sample_meta', type=str, default="/space/grp/rschwartz/rschwartz/cell_annotation_cortex.nf/mmus/2a/69d63b6c1090d6f10a71cc5662301c/GSE152715.1_sample_meta.tsv")
    if __name__ == "__main__":
        known_args, _ = parser.parse_known_args()
        return known_args

    
def main():
    # Parse command line arguments
    args = parse_arguments()
    # Set variables from arguments
    query_path = args.query_path
    assigned_celltypes_path = args.assigned_celltypes_path
    sample_meta = args.sample_meta
    markers_file = args.markers_file
    gene_mapping_path = args.gene_mapping 
    organism = args.organism
    
    gene_mapping = pd.read_csv(gene_mapping_path, sep=None, header=0)
    # Drop rows with missing values in the relevant columns
    gene_mapping = gene_mapping.dropna(subset=["ENSEMBL_ID", "OFFICIAL_SYMBOL"])
    # Set the index of gene_mapping to "ENSEMBL_ID" and ensure it's unique
    gene_mapping = gene_mapping.drop_duplicates(subset="ENSEMBL_ID")
    gene_mapping.set_index("ENSEMBL_ID", inplace=True) 

 

    # Load query and reference datasets
    study_name = os.path.basename(query_path).replace(".h5ad", "")
    assigned_celltypes = pd.read_csv(assigned_celltypes_path, sep=None, header=0)
    sample_meta = pd.read_csv(sample_meta, sep=None, header=0)
   # markers = pd.read_csv(markers_file, sep=None, header=0)
    os.makedirs(study_name, exist_ok=True)
    
    df = pd.read_csv(markers_file, sep=None, header=0)
    #write to mqc file
    df.to_csv(f"{study_name}/cell_type_markers_mqc.tsv", sep="\t", index=True)

    query = read_query(query_path, gene_mapping, new_meta=assigned_celltypes, sample_meta=sample_meta)
    query.obs.index = query.obs["index"]
    query.raw = query.copy()
    query = qc_preprocess(query)
    
    #plot_markers(query, markers_file, organism=organism)
    make_celltype_matrices(query, markers_file, organism=organism, study_name=study_name)
    
    
    query_subsets = {}
    for sample_name in query.obs["sample_name"].unique():
        query_subset = query[query.obs["sample_name"] == sample_name]
        
        query_subset = get_qc_metrics(query_subset, nmads=args.nmads)
        query_subsets[sample_name] = query_subset
        
        plot_umap_qc(query_subset, study_name=study_name, sample_name=sample_name)
        plot_jointplots(query_subset, study_name=study_name, sample_name=sample_name)
        
    # Count occurrences
    celltype_counts = (
        query.obs
        .groupby(["sample_name", "cell_type"])
        .size()                             # count cells per (sample, cell_type)
        .unstack(fill_value=0)              # pivot cell types into columns
        .reset_index()                      # make sample_name a column
    )
    celltype_counts.to_csv(os.path.join(study_name,"celltype_counts_mqc.tsv"), sep="\t", index=False)

    #combine query subsets
    query_combined = ad.concat(query_subsets.values(), axis=0) 
    ## make a table of counts by outliers
    outlier_counts = (
        query_combined.obs
        .groupby(["sample_name", "total_outlier"])
        .size()                             # count cells per (sample, cell_type)
        .unstack(fill_value=0)              # pivot cell types into columns
        .reset_index()                      # make sample_id a column
    )
    outlier_counts.to_csv(os.path.join(study_name,"outlier_counts_mqc.tsv"), sep="\t", index=False)

    
    
    
if __name__ == "__main__":
    main()
 
    
    

   