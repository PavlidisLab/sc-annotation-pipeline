#!/usr/bin/env python3

import warnings
warnings.filterwarnings("ignore")
import os
import random
import urllib.request
import pandas as pd
import numpy as np
import scanpy as sc
import anndata as ad
import cellxgene_census
import cellxgene_census.experimental


def setup(organism="homo_sapiens", version="2024-07-01"):
    organism = organism.replace(" ", "_")
    outdir = f"scvi-{organism}-{version}"

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    model_file_path = os.path.join(outdir, "model.pt")
    scvi_info = cellxgene_census.experimental.get_embedding_metadata_by_name(
        embedding_name="scvi",
        organism=organism,
        census_version=version,
    )

    model_link = scvi_info["model_link"]
    date = model_link.split("/")[5]
    url = os.path.join("https://cellxgene-contrib-public.s3.us-west-2.amazonaws.com/models/scvi/", date, organism, "model.pt")

    import ssl
    ctx = ssl._create_unverified_context()
    with urllib.request.urlopen(url, context=ctx) as response:
        with open(model_file_path, 'wb') as f:
            f.write(response.read())

    return outdir


def rename_cells(obs, rename_file="/space/grp/rschwartz/rschwartz/cell_annotation_cortex.nf/meta/rename_cells_mmus.tsv"):
    rename_df = pd.read_csv(rename_file, sep=None)
    rename_key = rename_df.columns[0]
    rename_cell_type = rename_df.columns[1]

    obs = obs[obs[rename_key].isin(rename_df[rename_key])]

    rename_mapping = dict(zip(rename_df[rename_key], rename_df[rename_cell_type]))
    ontology_mapping = dict(zip(rename_df[rename_cell_type], rename_df[rename_cell_type + "_uri"]))

    obs[rename_cell_type] = obs[rename_key].replace(rename_mapping)
    obs[rename_cell_type + "_uri"] = obs[rename_cell_type].map(ontology_mapping)

    return obs


def subsample_cells(data, filtered_ids, subsample=500, seed=42, organism="Homo sapiens",
                    rename_file="/space/grp/rschwartz/rschwartz/cell_annotation_cortex.nf/meta/rename_cells_mmus.tsv"):
    random.seed(seed)
    np.random.seed(seed)

    obs = data[data['soma_joinid'].isin(filtered_ids)]
    obs = rename_cells(obs, rename_file=rename_file)
    celltypes = obs["cell_type"].unique()
    print(obs["cell_type"].value_counts().reset_index())
    final_idx = []
    for celltype in celltypes:
        celltype_ids = obs[obs["cell_type"] == celltype]['soma_joinid'].values
        if len(celltype_ids) > subsample:
            subsampled_cell_idx = random.sample(list(celltype_ids), subsample)
        else:
            subsampled_cell_idx = celltype_ids.tolist()
        final_idx.extend(subsampled_cell_idx)

    return final_idx


def get_original_celltypes(columns_file, author_annotations_path):
    original_celltype_columns = pd.read_csv(columns_file, sep=None)

    original_celltypes = {}
    for file in os.listdir(author_annotations_path):
        if "obs.tsv" in file:
            dataset_title = file.split(".")[0]
            og_obs = pd.read_csv(os.path.join(author_annotations_path, file), sep=None)
            assert og_obs["observation_joinid"].nunique() == og_obs.shape[0]
            og_column = original_celltype_columns[original_celltype_columns["dataset_title"] == dataset_title]["author_cell_type"].values[0]
            og_obs["author_cell_type"] = og_obs[og_column]
            original_celltypes[dataset_title] = og_obs

    for dataset_title, obs in original_celltypes.items():
        original_celltypes[dataset_title]["new_observation_joinid"] = original_celltypes[dataset_title]["observation_joinid"].apply(lambda x: f"{dataset_title}_{x}")

    aggregate_obs = pd.concat([original_celltypes[ref_name] for ref_name in original_celltypes.keys()])
    duplicate_observation_joinid = aggregate_obs[aggregate_obs["new_observation_joinid"].duplicated()]
    duplicate_observation_joinid.columns
    assert aggregate_obs["new_observation_joinid"].nunique() == aggregate_obs.shape[0]

    return aggregate_obs


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
                 obs_filter=None, cell_columns=None, dataset_info=None, rename_file=None,
                 original_celltypes=None, seed=42):

    brain_cell_subsampled_ids = subsample_cells(cellxgene_obs_filtered, filtered_ids, subsample, seed=seed, organism=organism, rename_file=rename_file)
    adata = cellxgene_census.get_anndata(
        census=census,
        organism=organism,
        obs_value_filter=obs_filter,
        obs_column_names=cell_columns,
        obs_coords=brain_cell_subsampled_ids,
        var_value_filter="nnz > 10",
        obs_embeddings=["scvi"])

    sc.pp.filter_genes(adata, min_cells=3)
    sc.pp.filter_genes(adata, min_counts=200)

    print("Subsampling successful.")
    newmeta = adata.obs.merge(dataset_info, on="dataset_id", suffixes=(None, "y"))
    adata.obs = newmeta

    if isinstance(original_celltypes, pd.DataFrame) and not original_celltypes.empty:
        adata.obs = map_author_labels(adata.obs, original_celltypes)

    return adata


def get_cellxgene_obs(census, organism, organ="brain", primary_data=True, disease="normal"):
    value_filter = (
        f"tissue_general == '{organ}' and "
        f"is_primary_data == {str(primary_data)} and "
        f"disease == '{disease}'"
    )
    return cellxgene_census.get_obs(census, organism, value_filter=value_filter)


def get_census(census_version="2024-07-01", organism="homo_sapiens", subsample=5, assay=None, tissue=None, organ="brain",
               ref_collections=["Transcriptomic cytoarchitecture reveals principles of human neocortex organization"],
               original_celltypes=None,
               rename_file="/space/grp/rschwartz/rschwartz/cell_annotation_cortex.nf/meta/rename_cells.tsv", seed=42):

    census = cellxgene_census.open_soma(census_version=census_version)
    dataset_info = census.get("census_info").get("datasets").read().concat().to_pandas()

    cellxgene_obs = get_cellxgene_obs(census, organism, organ=organ, primary_data=True, disease="normal")

    cellxgene_obs = cellxgene_obs.merge(dataset_info, on="dataset_id", suffixes=(None, "_y"))
    cellxgene_obs.drop(columns=['soma_joinid_y'], inplace=True)
    cellxgene_obs_filtered = cellxgene_obs[cellxgene_obs['collection_name'].isin(ref_collections)]

    if assay:
        cellxgene_obs_filtered = cellxgene_obs_filtered[cellxgene_obs_filtered["assay"].isin(assay)]
    if tissue:
        cellxgene_obs_filtered = cellxgene_obs_filtered[cellxgene_obs_filtered["tissue"].isin(tissue)]

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

    filtered_ids = cellxgene_obs_filtered['soma_joinid'].values
    adata = extract_data(
        cellxgene_obs_filtered, filtered_ids,
        subsample=subsample, organism=organism,
        census=census, obs_filter=None,
        cell_columns=cell_columns, dataset_info=dataset_info, seed=seed,
        original_celltypes=original_celltypes,
        rename_file=rename_file
    )
    new_obs = rename_cells(adata.obs, rename_file=rename_file)
    new_adata = adata[adata.obs["soma_joinid"].isin(new_obs["soma_joinid"])].copy()
    new_obs_indexed = new_obs.set_index("soma_joinid")
    new_adata.obs = new_obs_indexed.loc[new_adata.obs["soma_joinid"].values].reset_index()

    for col in new_adata.obs.columns:
        if new_adata.obs[col].dtype.name == 'category':
            new_adata.obs[col] = pd.Categorical(new_adata.obs[col].cat.remove_unused_categories())

    return new_adata
