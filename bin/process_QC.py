#!/user/bin/python3
import warnings
warnings.filterwarnings("ignore")
from pathlib import Path
import os
import sys
import scanpy as sc
import numpy as np
import pandas as pd
import anndata as ad
import re
import warnings
from collections import defaultdict
import matplotlib.pyplot as plt
import seaborn as sns
import json
import argparse
import os
import json
from types import SimpleNamespace
from utils import *
from PIL import Image
import io
import os
import math
from upsetplot import UpSet, from_memberships


# Function to parse command line arguments
def parse_arguments():
    parser = argparse.ArgumentParser(description="Classify cells given 1 ref and 1 query")
    parser.add_argument('--organism', type=str, default='mus_musculus', help='Organism name (e.g., homo_sapiens)')
    parser.add_argument('--query_path', type=str, default="/space/grp/rschwartz/rschwartz/cell_annotation_cortex.nf/work/40/4adf027a41b7292db2847d7435c0f6/1373636_5M_Tim3_cKO.5XFAD_rep2_raw.h5ad")
    parser.add_argument('--assigned_celltypes_path', type=str, default="//space/grp/rschwartz/rschwartz/cell_annotation_cortex.nf/work/40/4adf027a41b7292db2847d7435c0f6/1373636_5M_Tim3_cKO.5XFAD_rep2_predicted_celltype.tsv")
    parser.add_argument('--markers_file', type=str, default="/space/grp/rschwartz/rschwartz/cell_annotation_cortex.nf/meta/cell_type_markers.tsv")
    parser.add_argument('--gene_mapping', type=str, default="/space/grp/rschwartz/rschwartz/cell_annotation_cortex.nf/meta/gemma_genes.tsv")
    parser.add_argument('--nmads',type=int, default=5)
    parser.add_argument('--sample_meta', type=str, default="/space/grp/rschwartz/rschwartz/cell_annotation_cortex.nf/work/40/4adf027a41b7292db2847d7435c0f6/GSE223423_sample_meta.tsv")
    parser.add_argument('--ref_keys', type=str, nargs="+", default=["subclass_cell_type","class_cell_type"], help="levels of granularity to classify corresponding to column names of rename_cells file")
    if __name__ == "__main__":
        known_args, _ = parser.parse_known_args()
        return known_args


def plot_joint_umap(query, study_name, sample_name):
    x_metric = "log1p_n_genes_by_counts"
    metrics = {
        "log1p_total_counts": ["counts_outlier", "umi_outlier", "genes_outlier"],
        "pct_counts_mito": "mito_outlier",
        "pct_counts_ribo": "ribo_outlier",
        "pct_counts_hb": "hb_outlier",
    }
    
    data = query.obs
    images = []
    for yval, hues in metrics.items():
        for hue in (hues if isinstance(hues, list) else [hues]):
            fig_joint = sns.jointplot(
                data=data,
                x=x_metric,
                y=yval,
                hue=hue,
                kind="scatter"
            )
        
        
            umap_fig = sc.pl.umap(
            query,
            color=hue,
            use_raw=False,
            save=None,
            show=False,
            title=f"{hue}",
            ncols=1,
            legend_loc="upper right",
            return_fig=True
            ) 


            joint_buf = io.BytesIO()
            fig_joint.savefig(joint_buf, format="png", bbox_inches='tight')
            plt.close(fig_joint.fig) 
            
            umap_buf = io.BytesIO()
            umap_fig.savefig(umap_buf, format="png", bbox_inches='tight')
            plt.close(umap_fig)
            
            joint_buf.seek(0)
            images.append(Image.open(joint_buf))
            
            umap_buf.seek(0)
            images.append(Image.open(umap_buf))
        
    scale = 0.5  # Resize to 50%
    resized_images = [img.resize((int(img.width * scale), int(img.height * scale))) for img in images]

    # Use resized dimensions
    img_width, img_height = resized_images[0].size
    grid_cols = 2
    grid_rows = math.ceil(len(resized_images) / grid_cols)

    combined_img = Image.new("RGB", (grid_cols * img_width, grid_rows * img_height), "white")

    for idx, img in enumerate(resized_images):
        row = idx // grid_cols
        col = idx % grid_cols
        x_offset = col * img_width
        y_offset = row * img_height
        combined_img.paste(img, (x_offset, y_offset))
    # replace slashes, spaces, weird stuff
    # fix this
    new_sample_name = re.sub(r'[^A-Za-z0-9._-]', '_', str(sample_name))
    out_path = os.path.join(study_name, f"{new_sample_name}_outliers_mqc.png")
    combined_img.save(out_path)




def plot_ct_umap(query, study_name):
    colors = [ref_keys[0],"leiden","sample_name"]
    
    fig = sc.pl.umap(
            query,
            color=colors,
            use_raw=False,
            show=False,
            title="",
            ncols=1,
            return_fig=True
        )

    out_path = os.path.join(study_name,"celltype_umap_mqc.png")
    fig.savefig(out_path, bbox_inches='tight')
    plt.close(fig)
 
def write_clc_files(query_combined, study_name, metrics=["counts_outlier", "mito_outlier", "ribo_outlier", "hb_outlier", "predicted_doublet", "umi_outlier", "genes_outlier"]):
    # change to long format
   # for metric in metrics:
    CLC_df = query_combined.obs[["sample_id", "cell_id"] + metrics].copy()

    
    CLC_df = CLC_df.melt(
        id_vars=["sample_id", "cell_id"],
        value_vars=metrics,
        value_name="value",
        var_name="category"
    )
        
    # change true and false to lower
    CLC_df["value"] = CLC_df["value"].astype(str).str.lower()
    # Save to TSV file
    CLC_df.to_csv(f"{study_name}_mask.tsv", sep="\t", index=False)


def plot_upset_by_group(obs, outlier_cols, group_col, outdir):
    os.makedirs(outdir, exist_ok=True)
    obs = obs.copy()
    obs["membership"] = obs[outlier_cols].apply(lambda row: tuple(c for c in outlier_cols if row[c]), axis=1)

    if group_col:
        for group in sorted(obs[group_col].unique()):
            counts = obs[obs[group_col] == group]["membership"].value_counts()
            data = from_memberships(counts.index, data=counts.values)
            plt.figure(figsize=(8, 4))
            UpSet(data, show_counts=True).plot()
            plt.suptitle(f"{group_col} = {group}")
            plt.tight_layout()
            new_group = re.sub(r'[^A-Za-z0-9._-]', '_', str(group))
            plt.savefig(os.path.join(outdir, f"{new_group}_upset_mqc.png"))
            plt.close()
            if counts.empty:
                continue
    else:   
        group_col = "study"
        group = "all"
        counts = obs["membership"].value_counts()
        if counts.empty:
            return
        data = from_memberships(counts.index, data=counts.values)
        plt.figure(figsize=(8, 4))
        UpSet(data, show_counts=True).plot()
        plt.suptitle(f"{group_col} = {group}")
        plt.tight_layout()
        plt.savefig(os.path.join(outdir, f"{group}_upset_mqc.png"))
        plt.close()
    

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
    ref_keys = args.ref_keys    
    
    # Load gene mapping file 
    gene_mapping = pd.read_csv(gene_mapping_path, sep=None, header=0)
    # Drop rows with missing values in the relevant columns
    gene_mapping = gene_mapping.dropna(subset=["ENSEMBL_ID", "OFFICIAL_SYMBOL"])
    # Set the index of gene_mapping to "ENSEMBL_ID" and ensure it's unique
    gene_mapping = gene_mapping.drop_duplicates(subset="ENSEMBL_ID")
    gene_mapping.set_index("ENSEMBL_ID", inplace=True) 

    # Load query and reference datasets
    study_name = os.path.basename(query_path).replace("_raw.h5ad", "")
    os.makedirs(study_name, exist_ok=True)

    assigned_celltypes = pd.read_csv(assigned_celltypes_path, sep=None, header=0)
    sample_meta = pd.read_csv(sample_meta, sep=None, header=0)
   # markers = pd.read_csv(markers_file, sep=None, header=0)
     
    query = read_query(query_path, gene_mapping, new_meta=assigned_celltypes, sample_meta=sample_meta)
    query.obs.index = query.obs["index"]
    query.raw = query.copy()
    make_celltype_matrices(query, markers_file, organism=organism, study_name=study_name, cell_type_key=assigned_celltypes.columns[2])

    query = qc_preprocess(query)
 
    query_subsets = {}
    for sample_name in query.obs["sample_name"].unique():
        query_subset = query[query.obs["sample_name"] == sample_name]
        
        query_subset = get_qc_metrics(query_subset, nmads=args.nmads)
        
        query_subsets[sample_name] = query_subset
        
       # plot_umap_qc(query_subset, study_name=study_name, sample_name=sample_name)
        plot_joint_umap(query_subset, study_name=study_name, sample_name=sample_name)
  
    #combine query subsets
    query_combined = ad.concat(query_subsets.values(), axis=0) 
     
    # Count occurrences
    celltype_counts = (
        query.obs
        .groupby(["sample_name", ref_keys[0]])
        .size()                             # count cells per (sample, cell_type)
        .unstack(fill_value=0)              # pivot cell types into columns
        .reset_index()                      # make sample_name a column
    )
    celltype_counts.to_csv(os.path.join(study_name,"celltype_counts_mqc.tsv"), sep="\t", index=False)
 
    plot_ct_umap(query_combined, study_name=study_name)
    
    cluster_celltypes = (
        query_combined.obs
        .groupby(["leiden", ref_keys[0]])
        .size() # count cells per (sample, cell_type)
        .unstack(fill_value=0)              # pivot cell types into columns
        .reset_index()                      # make sample_name a column
    )
    cluster_celltypes.to_csv(os.path.join(study_name,"cluster_celltypes_mqc.tsv"), sep="\t", index=False)

    # plot upset plots by sample and cell type
    
    outlier_cols = [
        "non_outlier",
        "counts_outlier", 
        "umi_outlier", 
        "genes_outlier",
        "mito_outlier", 
        "ribo_outlier", 
        "hb_outlier", 
        "predicted_doublet"
    ]
    plot_upset_by_group(query_combined.obs, outlier_cols, "sample_name", study_name)
    
    write_clc_files(query_combined, study_name, metrics=["counts_outlier", "mito_outlier", "ribo_outlier", "hb_outlier", "predicted_doublet", "umi_outlier", "genes_outlier"])
    
if __name__ == "__main__":
    main()
 
    
    

   