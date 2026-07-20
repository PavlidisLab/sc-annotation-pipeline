process DOWNLOAD_STUDIES {
    tag "$study_name"
    label 'process_single'

    input:
    val study_name
    val use_staging

    output:
    tuple val(study_name), path("${study_name}/"), emit: study_dir

    when:
    task.ext.when == null || task.ext.when

    script:
    def gemma_cmd = use_staging ? "gemma-cli-staging" : "gemma-cli"
    """
    ${gemma_cmd} getSingleCellDataMatrix \\
        -e ${study_name} \\
        --format mex \\
        --scale-type count \\
        --use-ensembl-ids \\
        -o ${study_name}
    """

    stub:
    """
    mkdir -p ${study_name}
    touch ${study_name}/matrix.mtx
    touch ${study_name}/barcodes.tsv
    touch ${study_name}/features.tsv
    """
}
