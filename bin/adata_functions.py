import warnings
warnings.filterwarnings("ignore")
import pandas as pd
import numpy as np
import scanpy as sc
import random
import cellxgene_census
import cellxgene_census.experimental
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


def setup(organism="homo_sapiens", version="2024-07-01"):
    organism=organism.replace(" ", "_") 
    #census = cellxgene_census.open_soma(census_version=version)
    outdir = f"scvi-{organism}-{version}"  # Concatenate strings using f-string
    
    # Create the output directory if it doesn't exist
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # Check if the model file exists
    model_file_path = os.path.join(outdir, "model.pt")
    #if not os.path.exists(model_file_path):
        # Get scVI model info
    scvi_info = cellxgene_census.experimental.get_embedding_metadata_by_name(
            embedding_name="scvi",
            organism=organism,
            census_version=version,
        )

        # Extract the model link
    model_link = scvi_info["model_link"]
    date = model_link.split("/")[5]
    url = os.path.join("https://cellxgene-contrib-public.s3.us-west-2.amazonaws.com/models/scvi/", date, organism, "model.pt")

    # Download the model using wget if it doesn't already exist
    subprocess.run(["wget", "--no-check-certificate", "-q", "-O", model_file_path, url])
# else:
     #   print(f"File already exists at {model_file_path}, skipping download.")

    return(outdir)

#clean up cellxgene ontologies

def rename_cells(obs, rename_file="/space/grp/rschwartz/rschwartz/cell_annotation_cortex.nf/meta/rename_cells_mmus.tsv"):
    # change to handle author_cell_type as first column name of rename df (distinguishes between author and cellxgene annotations)
    rename_df = pd.read_csv(rename_file, sep=None)
    rename_key = rename_df.columns[0]

    # Filter cells to keep only those with valid new cell types
    # this replaces the "restricted cell types" command line argument
    # not working for "glutamatergic neuron" in human
    
    obs = obs[obs[rename_key].isin(rename_df[rename_key])]
 
    # Create mapping dictionaries
    rename_mapping = dict(zip(rename_df[rename_key], rename_df['new_cell_type']))
    ontology_mapping = dict(zip(rename_df['new_cell_type'], rename_df['cell_type_ontology_term_id']))
    
    # Apply renaming
    obs['cell_type'] = obs[rename_key].replace(rename_mapping)  
    obs["cell_type_ontology_term_id"] = obs["cell_type"].map(ontology_mapping)
    
    return obs


# Subsample x cells from each cell type if there are n>x cells present
# calls "rename cells" to filter out generic/ambiguous cell types
# ensures equal representation of cell types in reference
def subsample_cells(data, filtered_ids, subsample=500, seed=42, organism="Homo sapiens", 
                    rename_file="/space/grp/rschwartz/rschwartz/cell_annotation_cortex.nf/meta/rename_cells_mmus.tsv"):
    """
    Subsample up to `subsample` cells per cell type, ensuring balanced representation.
    Filter and rename cells before sub-sampling (see rename_cells).
    """
    random.seed(seed)         # For `random`
    np.random.seed(seed)      # For `numpy`
    scvi.settings.seed = seed # For `scvi`
    
    # Filter data based on filtered_ids
    obs = data[data['soma_joinid'].isin(filtered_ids)]

    obs = rename_cells(obs, rename_file=rename_file)
    celltypes = obs["cell_type"].unique()
    print(obs["cell_type"].value_counts().reset_index())
    final_idx = []
    for celltype in celltypes:
        celltype_ids = obs[obs["cell_type"] == celltype]['soma_joinid'].values
        # Sample if there are enough observations, otherwise take all
        if len(celltype_ids) > subsample:
            subsampled_cell_idx = random.sample(list(celltype_ids), subsample)
        else:
            subsampled_cell_idx = celltype_ids.tolist()
        # Append subsampled indices to final list
        final_idx.extend(subsampled_cell_idx)

    # Return final indices
    return final_idx


def get_original_celltypes(columns_file="/space/grp/rschwartz/rschwartz/nextflow_eval_pipeline/meta/author_cell_annotations/original_celltype_columns.tsv", 
                           author_annotations_path="/space/grp/rschwartz/rschwartz/nextflow_eval_pipeline/meta/author_cell_annotations"):
    original_celltype_columns = pd.read_csv(columns_file, sep=None)

    original_celltypes = {}
    for file in os.listdir(author_annotations_path):
        if "obs.tsv" in file:
            dataset_title = file.split(".")[0]
            og_obs = pd.read_csv(os.path.join(author_annotations_path, file), sep=None)
            # check if all observation_joinid are unique
            assert og_obs["observation_joinid"].nunique() == og_obs.shape[0]
            og_column = original_celltype_columns[original_celltype_columns["dataset_title"] == dataset_title]["author_cell_type"].values[0]
            og_obs["author_cell_type"] = og_obs[og_column]
            original_celltypes[dataset_title] = og_obs

    for dataset_title, obs in original_celltypes.items():
        #original_celltypes[dataset_title]["new_dataset_title"] = dataset_title
        original_celltypes[dataset_title]["new_observation_joinid"] = original_celltypes[dataset_title]["observation_joinid"].apply(lambda x: f"{dataset_title}_{x}")
    
        # concat all original_celltypes
    aggregate_obs = pd.concat([original_celltypes[ref_name] for ref_name in original_celltypes.keys()])
    # find duplicate observation_joinid in aggregate_obs
    duplicate_observation_joinid = aggregate_obs[aggregate_obs["new_observation_joinid"].duplicated()]
    duplicate_observation_joinid.columns
    assert aggregate_obs["new_observation_joinid"].nunique() == aggregate_obs.shape[0]
   
    return aggregate_obs

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

def map_author_labels(obs, original_celltypes):
    obs["new_dataset_title"] = obs["dataset_title"].apply(lambda x: x.replace(" ", "_")
                                                                .replace("\\/", "_")
                                                                .replace("(", "")
                                                                .replace(")", "")
                                                                .replace("\\", "")
                                                                .replace("'", "")
                                                                .replace(":", "")
                                                                .replace(";", "")
                                                                .replace("&", ""))

    obs["new_observation_joinid"] = obs["new_dataset_title"].astype(str) + "_" + obs["observation_joinid"].astype(str)
    
    mapping = dict(zip(original_celltypes["new_observation_joinid"], original_celltypes["author_cell_type"]))
    obs["author_cell_type"] = obs["new_observation_joinid"].map(mapping)

    return obs

def extract_data(cellxgene_obs_filtered, filtered_ids, subsample=10, organism=None, census=None, 
    obs_filter=None, cell_columns=None, dataset_info=None, 
    original_celltypes=None, seed=42):
     
    brain_cell_subsampled_ids = subsample_cells(cellxgene_obs_filtered, filtered_ids, subsample, seed=seed, organism=organism)
    # Assuming get_seurat is defined to return an AnnData object
    adata = cellxgene_census.get_anndata(
        census=census,
        organism=organism,
        obs_value_filter=obs_filter,  # Ensure this is constructed correctly
        obs_column_names=cell_columns,
        obs_coords=brain_cell_subsampled_ids,
        var_value_filter = "nnz > 10",
        obs_embeddings=["scvi"])
    
    sc.pp.filter_genes(adata, min_cells=3)
    sc.pp.filter_genes(adata, min_counts=200)
    
    print("Subsampling successful.")
    newmeta = adata.obs.merge(dataset_info, on="dataset_id", suffixes=(None,"y"))
    adata.obs = newmeta
    
    if isinstance(original_celltypes, pd.DataFrame) and not original_celltypes.empty:
        adata.obs = map_author_labels(adata.obs, original_celltypes)
    # Assuming relabel_wrapper is defined
    # Convert all columns in adata.obs to factors (categorical type in pandas)
    return adata


def get_cellxgene_obs(census, organism, organ="brain", primary_data=True, disease="normal"):
    value_filter = (
        f"tissue_general == '{organ}' and "
        f"is_primary_data == {str(primary_data)} and "
        f"disease == '{disease}'"
    )
    return cellxgene_census.get_obs(census, organism, value_filter=value_filter)


def get_census(census_version="2024-07-01", organism="homo_sapiens", subsample=5, assay=None, tissue=None, organ="brain",
               ref_collections=["Transcriptomic cytoarchitecture reveals principles of human neocortex organization"," SEA-AD: Seattle Alzheimer’s Disease Brain Cell Atlas"], 
               original_celltypes = None,
               rename_file="/space/grp/rschwartz/rschwartz/cell_annotation_cortex.nf/meta/rename_cells.tsv",seed=42):

    census = cellxgene_census.open_soma(census_version=census_version)
    dataset_info = census.get("census_info").get("datasets").read().concat().to_pandas()
    
    cellxgene_obs = cellxgene_obs = get_cellxgene_obs(census, organism, organ=organ, primary_data=True, disease="normal")

    cellxgene_obs = cellxgene_obs.merge(dataset_info, on="dataset_id", suffixes=(None,"_y"))
    cellxgene_obs.drop(columns=['soma_joinid_y'], inplace=True)
    cellxgene_obs_filtered = cellxgene_obs[cellxgene_obs['collection_name'].isin(ref_collections)]
     
    if assay:
        cellxgene_obs_filtered = cellxgene_obs_filtered[cellxgene_obs_filtered["assay"].isin(assay)]
    if tissue:
        cellxgene_obs_filtered = cellxgene_obs_filtered[cellxgene_obs_filtered["tissue"].isin(tissue)]
         
    # Adjust organism naming for compatibility
    organism_name_mapping = {
        "homo_sapiens": "Homo sapiens",
        "mus_musculus": "Mus musculus"
    }
    organism = organism_name_mapping.get(organism, organism)

    cell_columns = [
        "assay", "cell_type", "cell_type_ontology_term_id", "tissue",
        "tissue_general", "suspension_type",
        "disease", "dataset_id", "development_stage",
        "soma_joinid", "observation_joinid"
    ]
    
    if isinstance(original_celltypes, pd.DataFrame) and not original_celltypes.empty:
       cellxgene_obs_filtered = map_author_labels(cellxgene_obs_filtered, original_celltypes)
        
    # Get embeddings for all data together
    filtered_ids = cellxgene_obs_filtered['soma_joinid'].values
    adata = extract_data(
        cellxgene_obs_filtered, filtered_ids,
        subsample=subsample, organism=organism,
        census=census, obs_filter=None,
        cell_columns=cell_columns, dataset_info=dataset_info, seed = seed,
        original_celltypes=original_celltypes
    )
    new_obs=rename_cells(adata.obs, rename_file=rename_file)
    new_adata = adata[adata.obs["soma_joinid"].isin(new_obs["soma_joinid"])].copy()
    new_obs = new_obs.set_index("soma_joinid")
    new_adata.obs = new_adata.obs.set_index("soma_joinid")
    new_adata.obs = new_adata.obs.loc[new_obs.index].copy()
    new_adata.obs = new_obs
     
   # for name, ref in refs.items(): 
    for col in new_adata.obs.columns:
        if new_adata.obs[col].dtype.name =='category':
            new_adata.obs[col] = pd.Categorical(new_adata.obs[col].cat.remove_unused_categories())
    
    return new_adata 



def process_query(query, model_file_path, batch_key="sample"):
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



def classify_cells(query, cutoff, probabilities):
    class_metrics = {}
    
    # Only use the first ref_key
    key = "cell_type" 
    #class_metrics[key] = {}

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
    
    # Store predictions and confidence in `query`
    query.obs[key] = predicted_classes
    #query["confidence"] = np.max(class_probs, axis=1)  # Store max probability as confidence
    
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
    outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
        np.median(M) + nmads * median_abs_deviation(M) < M
    )
    return outlier


def qc_preprocess(query):
    # check if any sample_id has fewer than 30 associated cwells
    sample_counts = query.obs["sample_id"].value_counts()
    if (sample_counts < 30).any():
        batch_key=None
    else:
        batch_key="sample_id"
    sc.pp.scrublet(query, batch_key=batch_key)
    # log normalize, comput neighbors and umap
    sc.pp.normalize_total(query, target_sum=1e4)
    sc.pp.log1p(query)
    sc.pp.highly_variable_genes(query, n_top_genes=2000, subset=False)
    sc.pp.pca(query)
    sc.pp.neighbors(query, n_neighbors=10, n_pcs=30)
    sc.tl.umap(query)
    sc.tl.leiden(query, resolution=0.3)
    
    return query



def mad(var, scale='normal'):
    """Median Absolute Deviation. Set scale='normal' for consistency with R's default."""
    med = np.median(var)
    mad = np.median(np.abs(var - med))
    if scale == 'normal':
        return mad * 1.4826  # for normally distributed data
    return mad


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
    query.var["mito"] = query.var["feature_name"].str.startswith(("MT", "mt", "Mt"))
    query.var["ribo"] = query.var["feature_name"].str.startswith(("RP", "Rp", "rp"))
    query.var["hb"] = query.var["feature_name"].str.startswith(("HB", "Hb","hb"))
    # fill NaN values with False
    query.var["mito"].fillna(False, inplace=True)
    query.var["ribo"].fillna(False, inplace=True)
    query.var["hb"].fillna(False, inplace=True) 

    sc.pp.calculate_qc_metrics(query, qc_vars=["mito", "ribo", "hb"], log1p=True, inplace=True, percent_top=[20], use_raw=True)

    metrics = {
        "log1p_total_counts": "umi_outlier",
        "log1p_n_genes_by_counts": "genes_outlier",
        "pct_counts_mito": "outlier_mito",
        "pct_counts_ribo": "outlier_ribo",
        "pct_counts_hb": "outlier_hb",
    }
    
    for metric, col_name in metrics.items():
        query.obs[col_name] = is_outlier(query, metric, nmads)

    lm_dict = get_lm(query, nmads=nmads)
    intercept = lm_dict["model"].params[0]
    slope = lm_dict["model"].params[1]
    

    query.obs["counts_outlier"] = (
        query.obs["log1p_n_genes_by_counts"] < (query.obs["log1p_total_counts"] * slope + (intercept - lm_dict["intercept_adjustment"]))
        ) | (
        query.obs["log1p_n_genes_by_counts"] > (query.obs["log1p_total_counts"] * slope + (intercept + lm_dict["intercept_adjustment"]))
        ) | (
        query.obs["umi_outlier"] ) | (query.obs["genes_outlier"])
        

    query.obs["total_outlier"] = (
        query.obs["counts_outlier"] | query.obs["outlier_mito"] | query.obs["outlier_ribo"] | query.obs["outlier_hb"] | query.obs["predicted_doublet"]
    )
    
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


def get_gene_to_celltype_map(markers_file, organism="mus_musculus"):
    # Read the marker file
    df = pd.read_csv(markers_file, sep="\t")
    df = df[df["markers"].notnull()]
    gene_to_celltype = {}

    for _, row in df.iterrows():
        cell_type = row["cell_type"]
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
        gene: f"{gene}: {'_'.join(set(celltypes))}"
        for gene, celltypes in gene_to_celltype.items()
    }
    return gene_ct_dict


def make_celltype_matrices(query, markers_file, organism="mus_musculus", study_name=""):
    # Drop vars with NaN feature names
    query = query[:, ~query.var["feature_name"].isnull()]
    query.var_names = query.var["feature_name"]
    
    #Make raw index match processed var index
    query.raw.var.index = query.raw.var["feature_name"]
    
    # Map cell type hierarchy
    query = map_celltype_hierarchy(query, markers_file=markers_file)

    # Read marker genes
    gene_ct_dict = get_gene_to_celltype_map(markers_file, organism=organism)
    # Collect all unique markers across all families/classes/cell types
    all_markers = list(gene_ct_dict.keys())
    valid_markers = [gene for gene in all_markers if gene in query.var_names]

    # Filter raw expression matrix to match query.var_names
    expr_matrix = query.raw.X.toarray()
    expr_matrix = pd.DataFrame(expr_matrix, index=query.obs.index, columns=query.raw.var.index)
    
    avg_expr = expr_matrix.groupby(query.obs["cell_type"]).mean()
    avg_expr = avg_expr.loc[:, valid_markers]
    
    # Scale expression across genes
    scaled_expr = (avg_expr - avg_expr.mean()) / avg_expr.std()
    scaled_expr = scaled_expr.loc[:, valid_markers]
    scaled_expr.fillna(0, inplace=True)

    # Rename columns: gene -> gene (celltype)
    scaled_expr.rename(columns=gene_ct_dict, inplace=True)

    # Save matrix
    os.makedirs(study_name, exist_ok=True)
    scaled_expr.to_csv(f"{study_name}/heatmap_mqc.tsv", sep="\t")

 