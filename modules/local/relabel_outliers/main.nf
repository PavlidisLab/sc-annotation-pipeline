process RELABEL_OUTLIERS {
    tag "$study_name"
    // label 'process_single'

    conda "/home/rschwartz/anaconda3/envs/scanpyenv"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://raschwaa/scanpyenv-minimal:latest' :
        'raschwaa/scanpyenv-minimal:latest' }"

    input:
    tuple val(study_name), path(celltype_file), path(mask_file)

    output:
    tuple val(study_name), path("*.relabeled.tsv"), emit: celltypes

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    python ${projectDir}/bin/relabel_outliers.py \\
        --celltype_file ${celltype_file} \\
        --mask_file ${mask_file} \\
        --output ${celltype_file.baseName}.relabeled.tsv
    """

    stub:
    """
    cp ${celltype_file} ${celltype_file.baseName}.relabeled.tsv
    """
}
