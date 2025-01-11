#!/user/bin/python3

from pathlib import Path
import os
import sys
import scanpy as sc
import numpy as np
import pandas as pd
import anndata as ad
import cellxgene_census
import scvi
from scipy.sparse import csr_matrix
import warnings
import cellxgene_census
import cellxgene_census.experimental
import scvi
from sklearn.ensemble import RandomForestClassifier
import adata_functions
from adata_functions import *
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
   # parser.add_argument('--organism', type=str, default='homo_sapiens', help='Organism name (e.g., homo_sapiens)')
  #  parser.add_argument('--census_version', type=str, default='2024-07-01', help='Census version (e.g., 2024-07-01)')
    parser.add_argument('--query_path', type=str, default="")
    parser.add_argument('--ref_path', type=str, default="/space/grp/rschwartz/rschwartz/nextflow_eval_pipeline/refs/whole_cortex.h5ad") #nargs ="+")
    parser.add_argument('--cutoff', type=float, default=0, help="Cutoff probability for classification, else cell will be assigned unknown")
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

    # Load query and reference datasets
    query = ad.read_h5ad(query_path)
    query_name = os.path.basename(query_path).replace(".h5ad", "")
    ref = ad.read_h5ad(ref_path, backed="r")

    # Fit a random forest classifier to the reference scvi embeddings and cell type annotations
    rfc = RandomForestClassifier(class_weight='balanced', random_state=SEED)
    rfc.fit(ref.obsm["scvi"], ref.obs["cell_type"].values)
    
    # Predict cell type using embeddings generated from scvi model
    probs = rfc.predict_proba(query.obsm["scvi"])
    prob_df = pd.DataFrame(probs, columns=rfc.classes_)
    query = classify_cells(query, cutoff, prob_df)
    mapping = dict(ref.obs[["cell_type", "cell_type_ontology_term_id"]].drop_duplicates().values)
    query.obs["cell_type_ontology_term_id"] = query.obs["cell_type"].map(mapping)
    
    os.makedirs(query_name, exist_ok=True)
    
    filtered_obs = query.obs[["sample_id","cell_id","cell_type", "cell_type_ontology_term_id"]]
    filtered_obs.to_csv(os.path.join(query_name,f"{query_name}_predicted_celltype.tsv"), sep="\t", index=False)

if __name__ == "__main__":
    main()
    
