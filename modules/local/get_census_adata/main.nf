process GET_CENSUS_ADATA {
    tag "$organism"
    label 'process_high'

    conda "/home/rschwartz/anaconda3/envs/scanpyenv"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'docker://raschwaa/census-minimal:latest' :
        'raschwaa/census-minimal:latest' }"

    input:
    val ref_collections
    val organism
    val organ
    val census_version
    val subsample_ref
    path rename_file
    val seed
    path original_celltype_columns
    path author_annotations_path

    output:
    path "refs/*.h5ad", emit: ref_paths

    when:
    task.ext.when == null || task.ext.when

    script:
    args = task.ext.args ?: ''
    original_cols_arg = original_celltype_columns ? "--original_celltype_columns ${original_celltype_columns}" : ''
    author_annot_arg = author_annotations_path ? "--author_annotations_path ${author_annotations_path}" : ''
    """
    python ${projectDir}/bin/get_census_adata.py \\
        --organism ${organism} \\
        --organ ${organ} \\
        --census_version ${census_version} \\
        --subsample ${subsample_ref} \\
        --ref_collections ${ref_collections} \\
        --rename_file ${rename_file} \\
        --seed ${seed} \\
        ${original_cols_arg} \\
        ${author_annot_arg} \\
        ${args}
    """

    stub:
    """
    mkdir -p refs
    touch refs/reference.h5ad
    """
}
