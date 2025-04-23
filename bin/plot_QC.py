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

# Function to parse command line arguments
def parse_arguments():
    parser = argparse.ArgumentParser(description="Classify cells given 1 ref and 1 query")
    parser.add_argument('--organism', type=str, default='mus_musculus', help='Organism name (e.g., homo_sapiens)')
  #  parser.add_argument('--census_version', type=str, default='2024-07-01', help='Census version (e.g., 2024-07-01)')
    parser.add_argument('--query_path', type=str, default="/space/grp/rschwartz/rschwartz/cell_annotation_cortex.nf/mmus/2a/69d63b6c1090d6f10a71cc5662301c/GSE152715.1.h5ad")
    parser.add_argument('--assigned_celltypes_path', type=str, default="/space/grp/rschwartz/rschwartz/cell_annotation_cortex.nf/mmus/2a/69d63b6c1090d6f10a71cc5662301c/GSE152715.1_predicted_celltype.tsv")
    parser.add_argument('--markers_table', type=str, default="/space/grp/rschwartz/rschwartz/cell_annotation_cortex.nf/meta/chatgpt_cell_type_markers.tsv")
    parser.add_argument('--rename_file', type=str, default = "/space/grp/Pipelines/sc-annotation-pipeline/meta/rename_cells_mmus.tsv")
    parser.add_argument('--nmads', type=int, default="5")
    parser.add_argument('--gene_mapping', type=str, default="/space/grp/rschwartz/rschwartz/cell_annotation_cortex.nf/meta/gemma_genes.tsv")
    parser.add_argument('--sample_meta', type=str, default="/space/grp/rschwartz/rschwartz/cell_annotation_cortex.nf/mmus/2a/69d63b6c1090d6f10a71cc5662301c/GSE152715.1_sample_meta.tsv")
    parser.add_argument('--gemma_username', type=str, default="raschwar")
    parser.add_argument('--gemma_password', type=str, default="7nddtt")
    if __name__ == "__main__":
        known_args, _ = parser.parse_known_args()
        return known_args


    
def read_query(query_path, gene_mapping, new_meta, sample_meta):
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
   
   
   
    sample_meta["sample_id"] = sample_meta["sample_id"].astype(str)
    query.obs = query.obs.merge(sample_meta, left_on="sample_id", right_on="sample_id", how="left", suffixes=("", "_y"))
    
    
    columns_to_drop = [col for col in query.obs.columns if col.endswith("_y")]
    query.obs.drop(columns=columns_to_drop, inplace=True)
    return query


def is_outlier(query, metric: str, nmads=3):
    M = query.obs[metric]
    outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
        np.median(M) + nmads * median_abs_deviation(M) < M
    )
    return outlier


def process_query(query):
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

    
def plot_jointplots(query, study_name, sample_name):
    os.makedirs(study_name, exist_ok=True)
    # Save query.obs to CSV
    query.obs.to_csv(f"{study_name}/{sample_name}_obs.tsv", sep="\t", index=False)
    tsv_path = os.path.abspath(f"{study_name}/{sample_name}_obs.tsv")

    # get path to Rscript
    rscript_path = os.path.join(os.path.dirname(__file__), "plot_jointplots.R")
    subprocess.run([
        "Rscript", rscript_path,
        tsv_path, study_name, sample_name
    ])

        

def plot_umap_qc(query, study_name, sample_name):
    colors = ["outlier_hb", "outlier_ribo", "outlier_mito","predicted_doublet","counts_outlier","total_outlier"]

    output_dir = os.path.join(study_name, sample_name)
    os.makedirs(output_dir, exist_ok=True)

    sc.pl.umap(
        query,
        color=colors,
        use_raw=False,
        save=None,
        show=False,
       # title=f"Sample {sample_name}",
        ncols=2)
        # Manually save the plot
    plt.savefig(os.path.join(output_dir, "umap_mqc.png"), dpi=150, bbox_inches='tight')
    plt.close()
            

def read_markers(markers_table, organism):
    df = pd.read_csv(markers_table, sep="\t", header=0)
    
    # Split markers column into list
    df['markers'] = df['markers'].str.split(',\s*', regex=True)

    # Build nested dict: family > class > cell_type
    nested_dict = defaultdict(lambda: defaultdict(dict))
    
    for _, row in df.iterrows():
        fam = row['family']
        cls = row['class']
        cell = row['cell_type']
        markers = row['markers']
        nested_dict[fam][cls][cell] = markers
    
    if organism == "mus_musculus":
        for fam in nested_dict:
            for cls in nested_dict[fam]:
                for cell in nested_dict[fam][cls]:
                    nested_dict[fam][cls][cell] = [
                        x.lower().capitalize() for x in nested_dict[fam][cls][cell]
                    ]

    return nested_dict

    
def map_celltype_hierarchy(query, markers_table):
    # Load the markers table
    df = pd.read_csv(markers_table, sep="\t", header=0)
    df.drop(columns="markers", inplace=True)
    query.obs = query.obs.merge(df, left_on="cell_type", right_on="cell_type", how="left", suffixes=("", "_y"))
    return query


def make_stable_colors(color_mapping_df):
    
    all_subclasses = sorted(color_mapping_df["new_cell_type"])
    # i need to hardcode a separate color palette based on the mmus mapping file
    # Generate unique colors for each subclass
    color_palette = sns.color_palette("husl", n_colors=len(all_subclasses))
    subclass_colors = dict(zip(all_subclasses, color_palette))
    return subclass_colors 


def plot_markers(query, markers_table, organism="mus_musculus"):
    # drop var with "nan" feature name
    query = query[:, ~query.var["feature_name"].isnull()]
    query.var_names = query.var["feature_name"]    
    query=map_celltype_hierarchy(query, markers_table=markers_table)
    nested_dict = read_markers(markers_table, organism)
    for family in query.obs["family"].unique():
        
        query_subset = query[query.obs["family"] == family]
        markers_to_plot = []
        for cls in nested_dict[family]:
            for cell_type in nested_dict[family][cls]:
                markers_to_plot.extend(nested_dict[family][cls][cell_type])
        # remove duplicates
        markers_to_plot = list(set(markers_to_plot))
        markers_to_plot = [marker for marker in markers_to_plot if marker in query.var_names]      
          
        if len(markers_to_plot) == 0:
            print(f"no markers to plot for {cell_type}")
            continue
        prefix = family.replace(" ", "_").replace("/", "_")
        plot_annotated_heatmap(query_subset, markers=markers_to_plot, groupby=["cell_type"], prefix=prefix)
        

def make_celltype_matrices(query, markers_table, organism="mus_musculus", study_name=""):
        # drop var with "nan" feature name
    query = query[:, ~query.var["feature_name"].isnull()]
    query.var_names = query.var["feature_name"]    
    query=map_celltype_hierarchy(query, markers_table=markers_table)
    nested_dict = read_markers(markers_table, organism)
    for family in query.obs["family"].unique():
        query_subset = query[query.obs["family"] == family].copy()
        markers_to_plot = []
        for cls in nested_dict[family]:
            for cell_type in nested_dict[family][cls]:
                markers_to_plot.extend(nested_dict[family][cls][cell_type])
        # remove duplicates
        markers_to_plot = list(set(markers_to_plot)) 
        prefix = family.replace(" ", "_").replace("/", "_")
            # sort query by groupby columns
        valid_markers = [gene for gene in markers_to_plot if gene in query_subset.var_names]
        # Extract expression data
        expr_matrix = query_subset.X.toarray()
        expr_matrix = pd.DataFrame(expr_matrix, index=query_subset.obs.index, columns=query_subset.var_names)
        # make a distance matrix for a heatmap from expr_matrix

        # Example: average expression per cell type
        avg_expr = expr_matrix.groupby(query_subset.obs['cell_type']).mean()
        scaled_expr = (avg_expr - avg_expr.mean()) / avg_expr.std()  
        scaled_expr = scaled_expr.loc[:, valid_markers]

        # Save to TSV for MultiQC heatmap input
        scaled_expr.to_csv(f"{study_name}/{prefix}_heatmap_mqc.tsv", sep='\t')

        
        
  
def plot_annotated_heatmap(query, markers, groupby=["cell_type"],
                           figsize=(5, 10), prefix=""):
    # sort query by groupby columns
    query = query[query.obs.sort_values(groupby).index] 
    valid_markers = [gene for gene in markers if gene in query.var_names]

    # Extract expression data
    expr_matrix =query[:, valid_markers].X.toarray()
    expr_matrix = pd.DataFrame(expr_matrix, index=query.obs.index, columns=valid_markers)
        # Extract categorical annotation
        
    annotations = query.obs[groupby].astype(str).copy()
    # Flatten unique categories from all groupby columns
    unique_categories = sorted(set(annotations.values.ravel()))

    # Generate colors
    num_colors = len(unique_categories)
    palette = sns.color_palette("husl", num_colors)  # Use a color palette with distinct colors 
    color_map = dict(zip(unique_categories, palette))

    # Create row_colors DataFrame
    row_colors = pd.DataFrame(index=annotations.index)

    legend_dict = {}
    for col in groupby:
        row_colors[col] = annotations[col].map(color_map)  # Process each column separately
        legend_dict[col] = [(label, color_map[label]) for label in sorted(annotations[col].unique())]

    # Map colors to annotations
    row_colors = row_colors.reset_index(drop=True)
    expr_matrix = expr_matrix.reset_index(drop=True)
    
    # Plot heatmap with annotations
    g = sns.clustermap(expr_matrix, 
                       row_colors=row_colors, 
                       row_cluster=False, 
                       col_cluster=False,
                       z_score=None,
                       standard_scale=1,
                       figsize=figsize)

    g.ax_heatmap.set_yticks([])
    g.ax_heatmap.set_yticklabels([])
    g.ax_heatmap.set_xlabel("")
    for label in g.ax_heatmap.get_xticklabels():
      label.set_rotation(90)
    g.ax_heatmap.set_title(prefix)
  #  plt.tight_layout()
    plt.savefig(f"{prefix}_heatmap.png", bbox_inches="tight")
        # Create a figure for the legends
    fig, axes = plt.subplots(len(legend_dict), 1, figsize=(6, len(legend_dict) * 1.5))
        # Loop through the legend dictionary and create the patches for each category
    for title, elements in legend_dict.items():
        patches = [Patch(facecolor=color, edgecolor="black", label=label) for label, color in elements]
        axes.legend(handles=patches, title=title, loc="upper left", frameon=True, ncol=3)
        axes.set_axis_off()

    plt.subplots_adjust(hspace=0.6)  
    plt.tight_layout()
    plt.savefig(f"{prefix}_legend.png")
        
def main():
    # Parse command line arguments
    args = parse_arguments()
    # Set variables from arguments
    query_path = args.query_path
    assigned_celltypes_path = args.assigned_celltypes_path
    sample_meta = args.sample_meta
    markers_table = args.markers_table
    gene_mapping_path = args.gene_mapping 
    organism = args.organism
    
    gene_mapping = pd.read_csv(gene_mapping_path, sep="\t", header=0)
    # Drop rows with missing values in the relevant columns
    gene_mapping = gene_mapping.dropna(subset=["ENSEMBL_ID", "OFFICIAL_SYMBOL"])
    # Set the index of gene_mapping to "ENSEMBL_ID" and ensure it's unique
    gene_mapping = gene_mapping.drop_duplicates(subset="ENSEMBL_ID")
    gene_mapping.set_index("ENSEMBL_ID", inplace=True)
    
   # subclass_colors = make_stable_colors(rename_cells_df)

    # Load query and reference datasets
    study_name = os.path.basename(query_path).replace(".h5ad", "")
    assigned_celltypes = pd.read_csv(assigned_celltypes_path, sep="\t", header=0)
    sample_meta = pd.read_csv(sample_meta, sep="\t", header=0)
   # markers = pd.read_csv(markers_table, sep="\t", header=0)
    os.makedirs(study_name, exist_ok=True)

    query = read_query(query_path, gene_mapping, new_meta=assigned_celltypes, sample_meta=sample_meta)
    query.obs.index = query.obs["index"]
    sc.pp.scrublet(query, batch_key="sample_id")
    query.raw = query.copy()
    query = process_query(query)
    
    #plot_markers(query, markers_table, organism=organism)
    make_celltype_matrices(query, markers_table, organism=organism, study_name=study_name)
    
    
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
 
    
    

   