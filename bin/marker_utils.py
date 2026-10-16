import warnings
warnings.filterwarnings("ignore")
import os
import pandas as pd
import seaborn as sns


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
