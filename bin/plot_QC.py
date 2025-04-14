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
    parser.add_argument('--query_path', type=str, default="/space/grp/rschwartz/rschwartz/cell_annotation_cortex.nf/work/9a/c8024bb30c691c937917bcee847f6b/GSE152715.1.h5ad")
    parser.add_argument('--assigned_celltypes_path', type=str, default="/space/grp/rschwartz/rschwartz/cell_annotation_cortex.nf/work/9a/c8024bb30c691c937917bcee847f6b/GSE152715.1_predicted_celltype.tsv")
    parser.add_argument('--markers_file', type=str, default="")
    parser.add_argument('--rename_file', type=str, default = "/space/grp/Pipelines/sc-annotation-pipeline/meta/rename_cells_mmus.tsv")
    parser.add_argument('--nmads', type=int, default="5")
    parser.add_argument('--gene_mapping', type=str, default="/space/grp/rschwartz/rschwartz/cell_annotation_cortex.nf/meta/gemma_genes.tsv")
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

def get_qc_metrics(query, nmads):
    query.var["mito"] = query.var["feature_name"].str.startswith(("MT", "mt", "Mt"))
    query.var["ribo"] = query.var["feature_name"].str.startswith(("RP", "Rp", "rp"))
    query.var["hb"] = query.var["feature_name"].str.startswith(("HB", "Hb","hb"))
    # fill NaN values with False
    query.var["mito"].fillna(False, inplace=True)
    query.var["ribo"].fillna(False, inplace=True)
    query.var["hb"].fillna(False, inplace=True) 

    sc.pp.calculate_qc_metrics(query, qc_vars=["mito", "ribo", "hb"], log1p=True, inplace=True, percent_top=[20], use_raw=False)

    metrics = {
        "log1p_total_counts": "outlier_total_counts",
        "log1p_n_genes_by_counts": "outlier_n_genes_by_counts",
        "pct_counts_mito": "outlier_mito",
        "pct_counts_ribo": "outlier_ribo",
        "pct_counts_hb": "outlier_hb",
        "pct_counts_in_top_20_genes": "outlier_top_20_genes",
    }
    
    for metric, col_name in metrics.items():
        query.obs[col_name] = is_outlier(query, metric, nmads)


    query.obs["counts_outlier"] = query.obs["outlier_total_counts"] | query.obs["outlier_n_genes_by_counts"]

    query.obs["total_outlier"] = (
        is_outlier(query, "log1p_total_counts", nmads) |
        is_outlier(query, "log1p_n_genes_by_counts", nmads) |
        is_outlier(query, "pct_counts_mito", nmads) |
        is_outlier(query, "pct_counts_ribo", nmads) |
        is_outlier(query, "pct_counts_hb", nmads) # |
        #is_outlier(query, "pct_counts_in_top_20_genes", nmads)
    )

    return query

    
def plot_jointplots(query, query_name):
    
    # First jointplot
    plot1 = sns.jointplot(
        data=query.obs,
        x="log1p_n_genes_by_counts",
        y="log1p_total_counts",
        hue="counts_outlier",
        kind="scatter"
    )
    # save jointplot to directory named after y val
    output_dir = os.path.join("multiqc_dir",query_name, "counts")
    os.makedirs(output_dir, exist_ok=True)
    plot1.savefig(os.path.join(output_dir, f"jointplot_mqc.png"))
    
    # Second jointplot
    plot2 = sns.jointplot(
        data=query.obs,
        x="log1p_n_genes_by_counts",
        y="pct_counts_mito",
        hue="outlier_mito",
        kind="scatter"
    )
    # save jointplot to directory named after y val
    output_dir = os.path.join("multiqc_dir",query_name, "mitochondrial")
    os.makedirs(output_dir, exist_ok=True)
    plot2.savefig(os.path.join(output_dir, "jointplot_mqc.png"))
    

    # Third jointplot
    plot3 = sns.jointplot(
        data=query.obs,
        x="log1p_n_genes_by_counts",
        y="pct_counts_ribo",
        hue="outlier_ribo",
        kind="scatter"
    )
    # save jointplot to directory named after y val
    output_dir = os.path.join("multiqc_dir",query_name, "ribosomal")
    os.makedirs(output_dir, exist_ok=True)
    plot3.savefig(os.path.join(output_dir, "jointplot_mqc.png"))
    
    
    plot4 = sns.jointplot(
        data=query.obs,
        x="log1p_n_genes_by_counts",
        y="pct_counts_hb",
        hue="outlier_hb",
        kind="scatter"
    )
    # save jointplot to directory named after y val
    output_dir = os.path.join("multiqc_dir",query_name, "hemoglobin")
    os.makedirs(output_dir, exist_ok=True)
    plot4.savefig(os.path.join(output_dir, "jointplot_mqc.png")) 
        
    

def plot_umap_qc(query, query_name, subclass_colors):
    colors = ["cell_type", "total_outlier", 
            "pct_counts_mito", "pct_counts_ribo", "pct_counts_hb","predicted_doublet"]

    dirname_mapping = {
        "cell_type": "celltype",
        "total_outlier": "total_outlier",
        "pct_counts_mito": "mitochondrial",
        "pct_counts_ribo": "ribosomal",
        "pct_counts_hb": "hemoglobin",
        "predicted_doublet": "predicted_doublet"
    }
    
    
    for color in colors:
        is_categorical = query.obs[color].dtype.name == "category" or query.obs[color].dtype == object
        dirname=dirname_mapping.get(color, None)
        # Set output directory
        output_dir = os.path.join("multiqc_dir",query_name, dirname)
        os.makedirs(output_dir, exist_ok=True)

        # Plot without saving
        sc.pl.umap(
            query,
            color=color,
            use_raw=False,
            save=None,
            show=False,
            title=f"{color} - {query_name}",
            palette=subclass_colors if is_categorical else None,
            color_map="viridis" if not is_categorical else None,
        )

        # Manually save the plot
        plt.savefig(os.path.join(output_dir, "umap_mqc.png"), dpi=150, bbox_inches='tight')
        plt.close()
            
        
def make_stable_colors(color_mapping_df):
    
    all_subclasses = sorted(color_mapping_df["new_cell_type"])
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
    gene_mapping_path = args.gene_mapping 
    rename_cells_df = pd.read_csv(args.rename_file, sep="\t", header=0)
    gene_mapping = pd.read_csv(gene_mapping_path, sep="\t", header=0)
    # Drop rows with missing values in the relevant columns
    gene_mapping = gene_mapping.dropna(subset=["ENSEMBL_ID", "OFFICIAL_SYMBOL"])
    # Set the index of gene_mapping to "ENSEMBL_ID" and ensure it's unique
    gene_mapping = gene_mapping.drop_duplicates(subset="ENSEMBL_ID")
    gene_mapping.set_index("ENSEMBL_ID", inplace=True)
    
    subclass_colors = make_stable_colors(rename_cells_df)

    # Load query and reference datasets
    query_name = os.path.basename(query_path).replace(".h5ad", "")
    
    assigned_celltypes = pd.read_csv(assigned_celltypes_path, sep="\t", header=0)
    
   # markers = pd.read_csv(markers_file, sep="\t", header=0)
    os.makedirs("multiqc_dir", exist_ok=True)

    query = read_adata(query_path, gene_mapping, assigned_celltypes)
    query.obs.index = query.obs["index"]
    sc.pp.scrublet(query, batch_key="sample_id")
    query.raw = query.copy()
    query = process_adata(query)

    for sample_id in query.obs["sample_id"].unique():
        query_subset = query[query.obs["sample_id"] == sample_id]
        
        query_subset = get_qc_metrics(query_subset, nmads=args.nmads)
        
        plot_umap_qc(query_subset, query_name=sample_id, subclass_colors=subclass_colors)
        plot_jointplots(query_subset, query_name=sample_id)
           
    # Count occurrences
    counts = (
        query.obs
        .groupby(["sample_id", "cell_type"])
        .size()                             # count cells per (sample, cell_type)
        .unstack(fill_value=0)              # pivot cell types into columns
        .reset_index()                      # make sample_id a column
    )
    counts.to_csv(os.path.join("multiqc_dir","celltype_counts_mqc.tsv"), sep="\t", index=False)

        # make a wide table for the query but with qc metrics
    
    # Add outlier column to query
    
    
if __name__ == "__main__":
    main()
 
    
    

   