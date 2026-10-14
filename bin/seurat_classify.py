#!/usr/bin/env python3
"""Turn Seurat label-transfer prediction scores into the per-level cell type tables written by scvi_classify.py.

Each cell gets its highest-scoring subclass (or "unknown" if that score is not above --cutoff), then labels are
aggregated to the higher levels with the same mapping file as the scVI classifiers."""
import argparse
import os
import warnings

import pandas as pd
import scanpy as sc

from classify_utils import classify_cells

warnings.filterwarnings("ignore")


def parse_arguments():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--query_path", required=True, help="processed query h5ad (cell order must match the scores)")
    p.add_argument("--scores", required=True, help="prediction_scores_seurat.tsv.gz from predict_seurat.R")
    p.add_argument("--cutoff", type=float, default=0, help="Cutoff on the top Seurat prediction score, else the cell is unknown")
    p.add_argument("--mapping_file", required=True)
    p.add_argument("--ref_keys", nargs="+", default=["subclass_cell_type", "class_cell_type", "family_cell_type"])
    known, _ = p.parse_known_args()
    return known


def main():
    args = parse_arguments()
    query_name = os.path.basename(args.query_path).replace(".h5ad", "")
    mapping_df = pd.read_csv(args.mapping_file, sep="\t")
    query = sc.read_h5ad(args.query_path, backed="r").obs.copy()
    prob_df = pd.read_csv(args.scores, sep="\t")
    if len(prob_df) != len(query):
        raise ValueError(f"{len(prob_df)} score rows but {len(query)} query cells")

    query = classify_cells(query, args.cutoff, prob_df, ref_keys=args.ref_keys, mapping_df=mapping_df)

    for key in args.ref_keys:
        mapping = dict(mapping_df[[key, f"{key}_uri"]].drop_duplicates().values)
        query[f"{key}_ontology_term_id"] = query[key].map(mapping)
        query[f"{key}_uri"] = "http://purl.obolibrary.org/obo/" + query[f"{key}_ontology_term_id"].str.replace(":", "_")
        query.drop(columns=[f"{key}_ontology_term_id"], inplace=True)
        out = query[["sample_id", "cell_id", key, f"{key}_uri"]].rename(columns={key: "cell_type", f"{key}_uri": "cell_type_uri"})
        out.to_csv(f"{query_name}_{key}.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
