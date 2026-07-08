import warnings
warnings.filterwarnings("ignore")
import pandas as pd
import numpy as np
import anndata as ad
import scvi
from sklearn.ensemble import RandomForestClassifier


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
