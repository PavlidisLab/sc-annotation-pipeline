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

import adata_functions
from adata_functions import *

import matplotlib.pyplot as plt
import seaborn as sns
import json
import argparse
import os
import json
from types import SimpleNamespace
from scipy.stats import median_abs_deviation


# Function to parse command line arguments
def parse_arguments():
    parser = argparse.ArgumentParser(description="Classify cells given 1 ref and 1 query")
   # parser.add_argument('--organism', type=str, default='homo_sapiens', help='Organism name (e.g., homo_sapiens)')
  #  parser.add_argument('--census_version', type=str, default='2024-07-01', help='Census version (e.g., 2024-07-01)')
    parser.add_argument('--query_path', type=str, default="/space/grp/rschwartz/rschwartz/cell_annotation_cortex.nf/work/ad/335c8108b0299526896436bdd33896/GSE199460.2.h5ad")
    parser.add_argument('--assigned_celltypes_path', type=str, default="/space/grp/rschwartz/rschwartz/cell_annotation_cortex.nf/results/mus_musculus_subsample_ref_500_2025-04-02_16-53-40/GSE199460.2/GSE199460.2_predicted_celltype.tsv")
    parser.add_argument('--markers_file', type=str, default="")
    if __name__ == "__main__":
        known_args, _ = parser.parse_known_args()
        return known_args


    
def read_adata(query_path, gene_mapping, new_meta):
    
    query = sc.read_h5ad(query_path)
    # filter query for cells and genes
    #sc.pp.filter_cells(query, min_counts=3)
    # sc.pp.filter_cells(query, min_genes =200)
    if "feature_name" not in query.var.columns:
        query.var = query.var.merge(gene_mapping["OFFICIAL_SYMBOL"], left_index=True, right_index=True, how="left")
        query.var.rename(columns={"OFFICIAL_SYMBOL": "feature_name"}, inplace=True)
        # make symbol the index
       # query.var.set_index("OFFICIAL_SYMBOL", inplace=True)
        #drop nan values
    else:
        query.var.set_index("feature_name", inplace=True)
    
    query.obs=query.obs.reset_index()
    query.obs["full_barcode"] = query.obs["sample_id"].astype(str) + "_" + query.obs["cell_id"].astype(str)
    new_meta["full_barcode"] = new_meta["sample_id"].astype(str) + "_" + new_meta["cell_id"].astype(str)
    
    query.obs = query.obs.merge(new_meta, left_on="full_barcode", right_on="full_barcode", how="left", suffixes=("", "_y"))
    columns_to_drop = [col for col in query.obs.columns if col.endswith("_y")]
    query.obs.drop(columns=columns_to_drop, inplace=True)
    return query


def is_outlier(adata, metric: str, nmads=3):
    M = adata.obs[metric]
    outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
        np.median(M) + nmads * median_abs_deviation(M) < M
    )
    return outlier


def process_adata(query):
    # log normalize, comput neighbors and umap
    sc.pp.normalize_total(query, target_sum=1e4)
    sc.pp.log1p(query)
    sc.pp.highly_variable_genes(query, n_top_genes=2000, subset=False)
    sc.pp.pca(query)
    sc.pp.neighbors(query, n_neighbors=10, n_pcs=30)
    sc.tl.umap(query)
    sc.tl.leiden(query)
    
    return query

def get_qc_metrics(query):
    query.var["mito"] = query.var["feature_name"].str.startswith(("MT", "mt", "Mt"))
    query.var["ribo"] = query.var["feature_name"].str.startswith(("RP", "Rp", "rp"))
    query.var["hb"] = query.var["feature_name"].str.startswith(("HB", "Hb","hb"))
    # fill NaN values with False
    query.var["mito"].fillna(False, inplace=True)
    query.var["ribo"].fillna(False, inplace=True)
    query.var["hb"].fillna(False, inplace=True) 

    sc.pp.calculate_qc_metrics(query, qc_vars=["mito", "ribo", "hb"], log1p=True, inplace=True)
    columns_to_describe = [
            col for col in query.obs.columns
            if re.search(r'pct|log1p|total_counts', col)
        ]

    df = query.obs[columns_to_describe].describe().T
    print(df)
    return query

    
def plot_jointplots(query, query_name):


    # First jointplot
    plot1 = sns.jointplot(
        data=query.obs,
        x="log1p_total_counts",
        y="log1p_n_genes_by_counts",
        hue="outlier",
        kind="scatter"
    )
    plot1.savefig(f"{query_name}_genes_vs_counts.png")

    # Second jointplot
    plot2 = sns.jointplot(
        data=query.obs,
        x="log1p_total_counts",
        y="log1p_total_counts_mito",
        hue="outlier",
        kind="scatter"
    )
    plot2.savefig(f"{query_name}_mito_counts.png")

    # Third jointplot
    plot3 = sns.jointplot(
        data=query.obs,
        x="log1p_total_counts",
        y="log1p_total_counts_ribo",
        hue="outlier",
        kind="scatter"
    )
    plot3.savefig(f"{query_name}_ribo_counts.png")
    
    plot4 = sns.jointplot(
        data=query.obs,
        x="log1p_total_counts",
        y="log1p_total_counts_hb",
        hue="outlier",
        kind="scatter"
    )
    plot4.savefig(f"{query_name}_hb_counts.png")
    
    

def plot_umap_qc(query, query_name, subclass_colors):
    colors = ["cell_type", "correct", "outlier", 
            "pct_counts_mito", "pct_counts_ribo", "pct_counts_hb", "log1p_n_genes_by_counts", "sample_id"]

    for color in colors:
        is_categorical = query.obs[color].dtype.name == "category" or query.obs[color].dtype == object

        sc.pl.umap(
            query,
            color=color,
            use_raw=False,
            save=f"_{query_name}_{color}_qc_umap.png",
            show=False,
            title=f"{color} - {query_name}",
            palette=subclass_colors if is_categorical else None,
            color_map="viridis" if not is_categorical else None,
        )

def make_stable_colors(color_mapping_df):
    
    all_subclasses = sorted(color_mapping_df["subclass"])
    # i need to hardcode a separate color palette based on the mmus mapping file
    # Generate unique colors for each subclass
    color_palette = sns.color_palette("husl", n_colors=len(all_subclasses))
    subclass_colors = dict(zip(all_subclasses, color_palette))
    return subclass_colors 
    
        
def main():
    SEED = 42
    random.seed(SEED)         # For `random`
    np.random.seed(SEED)      # For `numpy`
    # For `torch`'
    scvi.settings.seed = SEED # For `scvi`
    # Parse command line arguments
    args = parse_arguments()

    # Set variables from arguments
    query_path = args.query_path
    assigned_celltypes_path = args.assigned_celltypes_path
    markers_file = args.markers_file
    gene_mapping_path = "/space/grp/rschwartz/rschwartz/cell_annotation_cortex.nf/meta/gemma_genes.tsv" 

    gene_mapping = pd.read_csv(gene_mapping_path, sep="\t", header=0)
    # Drop rows with missing values in the relevant columns
    gene_mapping = gene_mapping.dropna(subset=["ENSEMBL_ID", "OFFICIAL_SYMBOL"])

    # Set the index of gene_mapping to "ENSEMBL_ID" and ensure it's unique
    gene_mapping = gene_mapping.drop_duplicates(subset="ENSEMBL_ID")
    gene_mapping.set_index("ENSEMBL_ID", inplace=True)
    # Load query and reference datasets
    query_name = os.path.basename(query_path).replace(".h5ad", "")
    assigned_celltypes = pd.read_csv(assigned_celltypes_path, sep="\t", header=0)
   # markers = pd.read_csv(markers_file, sep="\t", header=0)
    
    query = read_adata(query_path, gene_mapping, assigned_celltypes)

    query = process_adata(query)
    subclass_colors = make_stable_colors(assigned_celltypes)
    query = get_qc_metrics(query)
    plot_umap_qc(query, query_name, subclass_colors)
    plot_jointplots(query, query_name)
 
    # Add outlier column to query
    
    
    
 
    
    

   