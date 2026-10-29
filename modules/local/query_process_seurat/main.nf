process QUERY_PROCESS_SEURAT {
    tag "$query_name"
    label 'process_medium'

    conda "/home/rschwartz/anaconda3/envs/r4.3"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'docker://raschwaa/r4.3-seurat:latest' :
        'raschwaa/r4.3-seurat:latest' }"

    input:
    tuple val(study_name), val(query_name), path(h5ad_file)

    output:
    tuple val(study_name), val(query_name), path(h5ad_file), path("${h5ad_file.baseName}.rds"), emit: query_rds

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    Rscript ${projectDir}/bin/seurat_preprocessing.R \\
        --h5ad_file ${h5ad_file} \\
        --normalization_method ${params.seurat_normalization} \\
        --dims ${params.seurat_dims} \\
        --nfeatures ${params.seurat_nfeatures} \\
        ${args}
    """

    stub:
    """
    touch ${h5ad_file.baseName}.rds
    """
}
