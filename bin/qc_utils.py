import warnings
warnings.filterwarnings("ignore")
import numpy as np
import scanpy as sc
from scipy.stats import median_abs_deviation
from statsmodels.formula.api import ols


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
