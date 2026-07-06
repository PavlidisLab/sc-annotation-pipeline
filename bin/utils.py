import warnings
warnings.filterwarnings("ignore")
import pandas as pd
import numpy as np
import scanpy as sc
import random
import os
import anndata as ad
import scvi
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import *
from sklearn.preprocessing import label_binarize
from collections import defaultdict
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path
import subprocess
from scipy.stats import median_abs_deviation
from statsmodels.formula.api import ols
import warnings
# silence warnings
warnings.filterwarnings("ignore")


def relabel(adata, relabel_path, join_key="", sep=None):
    # Read the relabel table from the file
    relabel_df = pd.read_csv(relabel_path, sep=sep)  # Adjust the separator as needed
    # Take the first column as the join key
    if join_key=="":
        join_key = relabel_df.columns[0]
    # Ensure the join_key is in both the AnnData object and the relabel DataFrame
    if join_key not in adata.obs.columns:
        raise ValueError(f"{join_key} not found in AnnData object observations.")
    if join_key not in relabel_df.columns:
        raise ValueError(f"{join_key} not found in relabel DataFrame.")
    # Perform the left join to update the metadata
    adata.obs = adata.obs.merge(relabel_df, on=join_key, how='left', suffixes=(None, "_y"))
    columns_to_drop = [col for col in adata.obs.columns if col.endswith('_y')]
    adata.obs.drop(columns=columns_to_drop, inplace=True)
    return adata

def process_query(query, model_file_path, batch_key="sample", seed=42):
    scvi.settings.seed = seed  
    # Ensure the input AnnData object is valid
    if not isinstance(query, ad.AnnData):
        raise ValueError("Input must be an AnnData object.")

    # Assign ensembl_id to var
    #query.var["ensembl_id"] = query.var["feature_id"]
    if "feature_id" in query.var.columns:
        query.var.set_index("feature_id", inplace=True)
    query.obs["n_counts"] = query.X.sum(axis=1)
    query.obs["joinid"] = list(range(query.n_obs))
    query.obs["batch"] = query.obs[batch_key]

    # Filter out missing HGNC features
    #query = query[:, query.var["feature_name"].notnull().values].copy()

    # Prepare the query AnnData for scVI
    scvi.model.SCVI.prepare_query_anndata(query, model_file_path)
    vae_q = scvi.model.SCVI.load_query_data(query, model_file_path)

    # Set the model to trained and get latent representation
    vae_q.is_trained = True
    latent = vae_q.get_latent_representation()
    query.obsm["scvi"] = latent

    return query


def rfc_pred(ref, query, ref_keys, seed):
    """
    Fit a RandomForestClassifier at the most granular level and aggregate probabilities for higher levels.
    
    Parameters:
    - ref: Reference data with labels.
    - query: Query data for prediction.
    - ref_keys: List of ordered keys from most granular to highest level (e.g., ["rachel_subclass", "rachel_class", "rachel_family"]).
    - tree: Dictionary representing the hierarchy of classes.
    
    Returns:
    - probabilities: Dictionary with probabilities for each level of the hierarchy.
    """
    probabilities = {}
    
    # The most granular key is the first in the ordered list
    granular_key = ref_keys[0]
    
    # Initialize and fit the RandomForestClassifier at the most granular level
    rfc = RandomForestClassifier(class_weight='balanced', random_state=seed)
    rfc.fit(ref.obsm["scvi"], ref.obs[granular_key].values)
    # Predict probabilities at e most granular level
    probs_granular = rfc.predict_proba(query.obsm["scvi"])
    class_labels_granular = rfc.classes_
    base_score = rfc.score(query.obsm["scvi"], query.obs[granular_key].values)

    # Store granular level probabilities
    probabilities[granular_key] = {
        "probabilities": probs_granular,
        "class_labels": class_labels_granular,
        "accuracy": base_score
    }
    
    return probabilities



def classify_cells(query:pd.DataFrame, cutoff: float, probabilities: pd.DataFrame, ref_keys=["subclass_cell_type","class_cell_type"], mapping_df=None):
    
    # Only use the first ref_key
    # must be ordered from most granular to highest level
    key = ref_keys[0]

    # Extract the class labels and probabilities (DataFrame structure)
    class_labels = probabilities.columns.values  # Class labels are the column names
    class_probs = probabilities.values  # Probabilities as a numpy array
    
    predicted_classes = []
    
    if cutoff > 0:
        # Find the class with the maximum probability for each cell
        max_class_indices = np.argmax(class_probs, axis=1)  # Get the index of the max probability
        max_class_probs = np.max(class_probs, axis=1)  # Get the max probability
        
        # Set predicted classes to "unknown" if the max probability does not meet the threshold
        predicted_classes = [
            class_labels[i] if prob > cutoff else "unknown"
            for i, prob in zip(max_class_indices, max_class_probs)
        ]
    else:
        # Direct prediction without threshold filtering
        predicted_classes = class_labels[np.argmax(class_probs, axis=1)]
    
    # Store predictions in `query`
    query[key] = predicted_classes

    query = aggregate_labels(query=query, mapping_df=mapping_df, ref_keys=ref_keys, predicted=False)

    return query


def aggregate_labels(query: pd.DataFrame, mapping_df: pd.DataFrame, ref_keys: list, predicted=False):
    """
    Aggregate subclass labels or predicted labels to higher levels using mapping_df.
    If predicted=True, operates on columns like 'predicted_subclass', otherwise on obs columns.
    """
    for i in range(1, len(ref_keys)):
        lower_key = ref_keys[i-1]
        higher_level_key = ref_keys[i]
        mapping = mapping_df.set_index(lower_key)[higher_level_key].to_dict()
        if predicted:
            query["predicted_" + higher_level_key] = query["predicted_" + lower_key].map(mapping)
            query["predicted_" + higher_level_key] = query["predicted_" + higher_level_key].fillna(query["predicted_" + lower_key])
        else:
            query[higher_level_key] = query[lower_key].map(mapping)
            query[higher_level_key] = query[higher_level_key].fillna(query[lower_key])
    return query
    
    
# functions for QC plotting --------------------------

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
    median = np.median(M)
    mad = median_abs_deviation(M)

    upper_outlier = M > median + nmads * mad

    # For mito, only flag high values (low mito is fine)
    if metric == "pct_counts_mito":
        return upper_outlier

    lower_outlier = M < median - nmads * mad
    return upper_outlier | lower_outlier


def qc_preprocess(query):
    # add an option to compute doublets since this step is failing?
    
    # check if any sample_id has fewer than 30 associated cells
    sample_counts = query.obs["sample_id"].value_counts()
    if (sample_counts < 30).any():
        batch_key = None
    else:
        batch_key = "sample_id"
    try:
        sc.pp.scrublet(query, batch_key=batch_key)
    except Exception as e:
        print(f"scrublet failed: {e}")
        # add predicted_doublet column with all False
        query.obs["predicted_doublet"] = False
    # log normalize, comput neighbors and umap
    sc.pp.normalize_total(query, target_sum=1e4)
    sc.pp.log1p(query)
    sc.pp.highly_variable_genes(query, n_top_genes=2000, subset=False)
    sc.pp.pca(query)
    sc.pp.neighbors(query, n_neighbors=10, n_pcs=30)
    sc.tl.umap(query)
    sc.tl.leiden(query, resolution=0.3)
    
    return query


def get_lm(query, nmads=5, scale="normal"):
    # Assume dataset is an AnnData object
    # Fit linear model: log10(n_genes_per_cell) ~ log10(counts_per_cell)
    lm_model = ols(formula='log1p_n_genes_by_counts ~ log1p_total_counts', data=query.obs).fit()
    # Calculate residuals
    residuals = lm_model.resid
    # If data is normally distributed, this is similar to std 
    mad_residuals = median_abs_deviation(residuals, scale=scale)
    # Intercept adjustment (add for upper bound, subtract for lower bound)
    intercept_adjustment = np.median(residuals) + nmads * mad_residuals
    return {
        "model": lm_model,
        "intercept_adjustment": intercept_adjustment
    }
    
    
def get_qc_metrics(query, nmads):
    """
    Calculate QC metrics and flag outliers.

    Parameters:
        query: AnnData object
        nmads: dict with keys 'mito', 'umi', 'genes', 'counts' for per-metric NMAD values.
    """

    query.var["mito"] = query.var["feature_name"].str.startswith(("MT", "mt", "Mt"))
    query.var["ribo"] = query.var["feature_name"].str.startswith(("RP", "Rp", "rp"))
    query.var["hb"] = query.var["feature_name"].str.startswith(("HB", "Hb","hb"))
    # fill NaN values with False
    query.var["mito"].fillna(False, inplace=True)
    query.var["ribo"].fillna(False, inplace=True)
    query.var["hb"].fillna(False, inplace=True)

    sc.pp.calculate_qc_metrics(query, qc_vars=["mito", "ribo", "hb"], log1p=True, inplace=True, percent_top=[20], use_raw=True)

    metrics = {
        "umi": "log1p_total_counts",
        "genes": "log1p_n_genes_by_counts",
        "mito": "pct_counts_mito",
    }

    for key, metric in metrics.items():
        query.obs[f"{key}_outlier"] = is_outlier(query, metric, nmads[key])

    lm_dict = get_lm(query, nmads=nmads["counts"])
    intercept = lm_dict["model"].params[0]
    slope = lm_dict["model"].params[1]
    

    query.obs["counts_outlier"] = (
        query.obs["log1p_n_genes_by_counts"] < (query.obs["log1p_total_counts"] * slope + (intercept - lm_dict["intercept_adjustment"]))
        ) | (
        query.obs["log1p_n_genes_by_counts"] > (query.obs["log1p_total_counts"] * slope + (intercept + lm_dict["intercept_adjustment"]))
        ) #| (
      #  query.obs["umi_outlier"] ) | (query.obs["genes_outlier"])
        

    metrics = ["counts_outlier", "mito_outlier", "ribo_outlier", "hb_outlier", "umi_outlier", "genes_outlier", "predicted_doublet"]
    # Ensure all metrics exist as boolean columns
    existing_metrics = [m for m in metrics if m in query.obs.columns] 
    # Calculate total_outlier efficiently
    query.obs["total_outlier"] = query.obs[existing_metrics].any(axis=1)
    
    query.obs["non_outlier"] = ~query.obs["total_outlier"]

    return query

            
def map_celltype_hierarchy(query, markers_file):
    # Load the markers table
    df = pd.read_csv(markers_file, sep=None, header=0)
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


def get_gene_to_celltype_map(df, organism="mus_musculus"):
    # Read the marker file
    df = df[df["markers"].notnull()]
    gene_to_celltype = {}

    for _, row in df.iterrows():
        cell_type = row["shortname"]
        genes = [gene.strip() for gene in row["markers"].split(",")]
        for gene in genes:
            if organism == "mus_musculus":
                gene = gene.lower().capitalize()
                # Handle multiple cell types mapping to the same gene
            if gene not in gene_to_celltype:
                gene_to_celltype[gene] = []
            gene_to_celltype[gene].append(cell_type)

    # Join multiple cell types into one label if needed
    gene_ct_dict = {
        gene: f"{'_'.join(set(celltypes))}: {gene}"
        for gene, celltypes in gene_to_celltype.items()
    }
    return gene_ct_dict


def make_celltype_matrices(query, markers_file, organism="mus_musculus", study_name="", cell_type_key="subclass_cell_type"):
    if cell_type_key not in query.obs.columns:
        return

    # Drop vars with NaN feature names
    query = query[:, ~query.var["feature_name"].isnull()]
    query.var_names = query.var["feature_name"]

    markers_df = pd.read_csv(markers_file, sep="\t")
    markers_df = markers_df[markers_df["organism"] == organism]
    level = cell_type_key.replace("_cell_type", "")
    if "level" in markers_df.columns:
        markers_df = markers_df[markers_df["level"] == level]
    ontology_mapping = markers_df.set_index("cell_type")["shortname"].to_dict()

    query.raw.var.index = query.raw.var["feature_name"]

    gene_ct_dict = get_gene_to_celltype_map(markers_df, organism=organism)
    all_markers = list(gene_ct_dict.keys())
    valid_markers = [gene for gene in all_markers if gene in query.var_names]

    expr_matrix = query.raw.X.toarray()
    expr_matrix = pd.DataFrame(expr_matrix, index=query.obs.index, columns=query.raw.var.index)

    avg_expr = expr_matrix.groupby(query.obs[cell_type_key]).mean()
    avg_expr = avg_expr.loc[:, valid_markers]

    scaled_expr = (avg_expr - avg_expr.mean()) / avg_expr.std()
    scaled_expr = scaled_expr.loc[:, valid_markers]
    scaled_expr.fillna(0, inplace=True)

    scaled_expr.rename(columns=gene_ct_dict, inplace=True)
    sorted_columns = sorted(scaled_expr.columns, key=lambda x: x.split(":")[0])
    scaled_expr = scaled_expr[sorted_columns]

    overlap = list(set(markers_df["cell_type"]).intersection(scaled_expr.index))
    sorted_cell_types = sorted(
        overlap,
        key=lambda x: ontology_mapping[x] if not pd.isna(ontology_mapping.get(x)) else x
    )

    scaled_expr = scaled_expr.loc[sorted_cell_types, :]

    os.makedirs(study_name, exist_ok=True)
    scaled_expr.to_csv(os.path.join(study_name, f"{study_name}_{cell_type_key}_heatmap_mqc.tsv"), sep="\t")

 