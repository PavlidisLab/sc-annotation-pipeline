"""Mito genes are found by Ensembl id on Gemma queries, whose chrM symbols are bare names (COX1, ND4).

Run: /home/rschwartz/anaconda3/envs/scanpyenv/bin/python -m pytest -s tests/test_mito_genes.py
Uses one human and one mouse Gemma query sample from get_gemma_data.nf (local files, no Gemma calls);
tests skip if the files are not on disk. Override paths with MITO_TEST_HS / MITO_TEST_MM.
"""
import os
import sys

import numpy as np
import pandas as pd
import pytest

BIN = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "bin")
sys.path.insert(0, BIN)
import scanpy as sc  # noqa: E402
from qc_utils import MITO_GENE_FILES, flag_mito_genes, get_qc_metrics, read_query  # noqa: E402

GEMMA = "/space/grp/rschwartz/rschwartz/get_gemma_data.nf"
QUERIES = {
    "hs": os.environ.get("MITO_TEST_HS", f"{GEMMA}/all_homo_sapiens_samples/h5ad/GSE157827/GSE157827_527211_GSM4775579.h5ad"),
    "mm": os.environ.get("MITO_TEST_MM", f"{GEMMA}/study_names_mouse.txt_author_true_process_samples_true/h5ad/GSE185454/GSE185454_1010394_GSM5615602.h5ad"),
}
NMADS = {"mito": 20, "umi": 5, "genes": 5, "counts": 5}


def chrm(org):
    return pd.read_csv([p for p in MITO_GENE_FILES if p.endswith(f"_{org}.tsv")][0], sep="\t", dtype=str)


@pytest.mark.parametrize("org", ["hs", "mm"])
def test_mito_found_on_gemma_query(org):
    path = QUERIES[org]
    if not os.path.exists(path):
        pytest.skip(f"{path} not on disk")
    raw = sc.read_h5ad(path)
    meta = raw.obs[["sample_id", "cell_id"]].astype(str).copy()
    sample_meta = pd.DataFrame({"sample_id": meta["sample_id"].unique(), "sample_name": meta["sample_id"].unique()})
    q = read_query(path, gene_mapping=None, new_meta=meta, sample_meta=sample_meta)  # feature_name -> var index
    q.raw = q.copy()

    ref = chrm(org)
    # the bug: the old name rule finds no chrM gene on Gemma symbols
    old = q.var["feature_name"].astype(str).str.startswith(("MT", "mt", "Mt"))
    assert not q.var.loc[old, "ensembl_id"].isin(ref["ensembl_id"]).any()

    mask, rule = flag_mito_genes(q.var)
    assert rule == "ensembl_id"
    found = q.var.loc[mask.values, "ensembl_id"]
    assert found.isin(ref["ensembl_id"]).all()          # no nuclear genes
    assert len(found) >= 13                              # the 13 protein-coding chrM genes
    q = get_qc_metrics(q, nmads=NMADS)
    pct = q.obs["pct_counts_mito"]
    print(f"{org}: {os.path.basename(path)} rule={rule} n_mito_genes={len(found)} "
          f"old_rule_genes={int(old.sum())} median_pct_mito={np.median(pct):.2f}")
    assert np.median(pct) > 0


def test_name_fallback_without_ensembl_ids():
    var = pd.DataFrame({"feature_name": ["MT-CO1", "mt-Nd4", "MTOR", "GAPDH"]}, index=["MT-CO1", "mt-Nd4", "MTOR", "GAPDH"])
    mask, rule = flag_mito_genes(var)
    assert rule.startswith("feature_name")
    assert mask.tolist() == [True, True, False, False]
