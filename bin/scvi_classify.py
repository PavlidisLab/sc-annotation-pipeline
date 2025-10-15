#!/user/bin/python3

from pathlib import Path
import os
import sys
import scanpy as sc
import numpy as np
import pandas as pd
import anndata as ad
import cellxgene_census
from scipy.sparse import csr_matrix
import warnings
import cellxgene_census
import cellxgene_census.experimental
import scvi
from sklearn.ensemble import RandomForestClassifier
from utils import *
from pathlib import Path
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, classification_report, confusion_matrix, roc_curve, auc
from sklearn.preprocessing import label_binarize
import matplotlib.pyplot as plt
import seaborn as sns
import json
import argparse
import os
import json
from types import SimpleNamespace

# Function to parse command line arguments
def parse_arguments():
    parser = argparse.ArgumentParser(description="Classify cells given 1 ref and 1 query")
    parser.add_argument('--query_path', type=str, default="")
    parser.add_argument('--ref_path', type=str, default="/space/grp/rschwartz/rschwartz/nextflow_eval_pipeline/refs/whole_cortex.h5ad") #nargs ="+")
    parser.add_argument('--cutoff', type=float, default=0, help="Cutoff probability for classification, else cell will be assigned unknown")
    parser.add_argument('--ref_keys', type=str, nargs="+", default=["subclass_cell_type","class_cell_type"], help="levels of granularity to classify corresponding to column names of rename_cells file")
    parser.add_argument('--mapping_df', type=str, default="/space/grp/Pipelines/sc-annotation-pipeline/meta/rename_cells_mmus_author.tsv", help="cell type taxonomy mapping file")
    if __name__ == "__main__":
        known_args, _ = parser.parse_known_args()
        return known_args

    
    
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
    ref_path = args.ref_path
    cutoff = args.cutoff
    ref_keys = args.ref_keys
    mapping_df = pd.read_csv(args.mapping_df, sep='\t')

    # Load query and reference datasets
    query_h5ad = sc.read_h5ad(query_path)
    query_name = os.path.basename(query_path).replace(".h5ad", "")
    ref = ad.read_h5ad(ref_path, backed="r")

    # Fit a random forest classifier to the reference scvi embeddings and cell type annotations
    # Training on the subclass level of granularity
    rfc = RandomForestClassifier(class_weight='balanced', random_state=SEED)
    rfc.fit(ref.obsm["scvi"], ref.obs["subclass_cell_type"].values)
    
    # Predict cell type using embeddings generated from scvi model
    probs = rfc.predict_proba(query_h5ad.obsm["scvi"])
    prob_df = pd.DataFrame(probs, columns=rfc.classes_)


    query = query_h5ad.obs
    # Classify cells based on probability cutoff and aggregate results at different levels of granularity
    query = classify_cells(query, cutoff, prob_df, ref_keys=ref_keys, mapping_df=mapping_df)
    
    # map to ontology terms at each level
    for key in ref_keys:
        mapping = dict(mapping_df[[key, f"{key}_uri"]].drop_duplicates().values)
        query[f"{key}_ontology_term_id"] = query[key].map(mapping)
        query[f"{key}_uri"] = f"http://purl.obolibrary.org/obo/" + query[f"{key}_ontology_term_id"].str.replace(":","_")
        # drop original ontology term id column
        query.drop(columns=[f"{key}_ontology_term_id"], inplace=True)
        
    os.makedirs(query_name, exist_ok=True)

    columns_to_keep = ["sample_id", "cell_id"] + [key for key in ref_keys] + [f"{key}_uri" for key in ref_keys]
    filtered_obs = query[columns_to_keep]
    filtered_obs.to_csv(f"{query_name}_predicted_celltype.tsv", sep="\t", index=False)

if __name__ == "__main__":
    main()
    
