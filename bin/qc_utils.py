import warnings
warnings.filterwarnings("ignore")
import os
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import median_abs_deviation

# chrM genes per organism (Ensembl ID, nomenclature symbol, NCBI ID, Gemma symbol), shipped with the pipeline
MITO_GENE_FILES = [
    os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "assets", f"mito_genes_{org}.tsv")
    for org in ("hs", "mm")
]
ENSEMBL_GENE_RE = r"^ENS[A-Z]*G\d{11}$"
ENSEMBL_ID_COLUMNS = ("ensembl_id", "gene_id", "gene_ids", "feature_id")


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
        # fall back to the Ensembl id for genes with no symbol, so var_names stays
        # all-string (a NaN float breaks var_names_make_unique's string concatenation).
        # astype(object) first: feature_name is often a pandas Categorical, and fillna
        # can't insert values that aren't already in the category list.
        query.var["feature_name"] = query.var["feature_name"].astype(object).fillna(query.var.index.to_series())
        # keep the original (Ensembl) ids: get_qc_metrics flags mito genes by id, since Gemma symbols
        # for chrM genes are bare names (COX1, ND4) that no MT- prefix rule can find
        if "ensembl_id" not in query.var.columns:
            query.var["ensembl_id"] = query.var.index.astype(str)
        query.var.set_index("feature_name", inplace=True, drop=False)
        query.var_names_make_unique()

    query.obs=query.obs.reset_index()
    query.obs["full_barcode"] = query.obs["sample_id"].astype(str) + "_" + query.obs["cell_id"].astype(str)
    new_meta["full_barcode"] = new_meta["sample_id"].astype(str) + "_" + new_meta["cell_id"].astype(str)

    query.obs = query.obs.merge(new_meta, left_on="full_barcode", right_on="full_barcode", how="left", suffixes=("", "_y"))

    sample_meta["sample_id"] = sample_meta["sample_id"].astype(str)
    query.obs = query.obs.merge(sample_meta, left_on="sample_id", right_on="sample_id", how="left", suffixes=("", "_y"))

    # use_gemma=false means no metadata to join, so sample_name is all NaN - fall back to sample_id.
    # astype(object) first: sample_name is often a pandas Categorical (h5ad-sourced data
    # commonly stores string columns this way), and fillna can't insert values that aren't
    # already in the category list.
    query.obs["sample_name"] = query.obs["sample_name"].astype(object).fillna(query.obs["sample_id"])

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
    x = query.obs["log1p_total_counts"].to_numpy()
    y = query.obs["log1p_n_genes_by_counts"].to_numpy()
    slope, intercept = np.polyfit(x, y, 1)
    # Calculate residuals
    residuals = y - (intercept + slope * x)
    # If data is normally distributed, this is similar to std
    mad_residuals = median_abs_deviation(residuals, scale=scale)
    # Intercept adjustment (add for upper bound, subtract for lower bound)
    intercept_adjustment = np.median(residuals) + nmads * mad_residuals
    return {
        "slope": slope,
        "intercept": intercept,
        "intercept_adjustment": intercept_adjustment
    }


def get_ensembl_ids(var):
    """Ensembl gene ids aligned to var (version suffix stripped), from var_names or an id column; None if absent."""
    candidates = [pd.Series(var.index, index=var.index)]
    candidates += [var[c] for c in ENSEMBL_ID_COLUMNS if c in var.columns]
    for ids in candidates:
        ids = ids.astype(str).str.replace(r"\.\d+$", "", regex=True)
        if ids.str.match(ENSEMBL_GENE_RE).mean() > 0.5:
            return ids
    return None


def load_mito_ids(paths=MITO_GENE_FILES):
    return set(pd.concat([pd.read_csv(p, sep="\t", dtype=str) for p in paths])["ensembl_id"])


def flag_mito_genes(var, mito_ids=None):
    """
    Boolean mask of mitochondrially encoded (chrM) genes in var, and the rule used.

    Matches Ensembl ids against assets/mito_genes_{hs,mm}.tsv. NCBI's main symbol for these genes is the
    bare name (COX1, ND4), which Gemma imports, so a name prefix rule finds no chrM gene on Gemma queries
    and many nuclear genes (MTOR, MTHFR, ...). The name rule (MT-/mt- prefix) is used only when var has
    no Ensembl ids.
    """
    ids = get_ensembl_ids(var)
    if ids is not None:
        mask = ids.isin(load_mito_ids() if mito_ids is None else mito_ids).values
        rule = "ensembl_id"
    else:
        mask = var["feature_name"].astype(str).str.match(r"^(MT|mt)-").values
        rule = "feature_name prefix MT-/mt- (no Ensembl ids in var)"
    print(f"mito genes: rule={rule}, n_matched={int(mask.sum())}")
    if mask.sum() == 0:
        print("WARNING: no mitochondrial genes matched; pct_counts_mito will be 0 for every cell")
    return pd.Series(mask, index=var.index), rule


def get_qc_metrics(query, nmads):
    """
    Calculate QC metrics and flag outliers.

    Parameters:
        query: AnnData object
        nmads: dict with keys 'mito', 'umi', 'genes', 'counts' for per-metric NMAD values.
    """

    query.var["mito"], query.uns["mito_rule"] = flag_mito_genes(query.var)
    query.var["ribo"] = query.var["feature_name"].str.startswith(("RP", "Rp", "rp"))
    query.var["hb"] = query.var["feature_name"].str.startswith(("HB", "Hb","hb"))
    # fill NaN values with False
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
    intercept = lm_dict["intercept"]
    slope = lm_dict["slope"]


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
