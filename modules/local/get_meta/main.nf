process GET_META {
    tag "$study_name"
    // label 'process_single'

    conda "/home/rschwartz/anaconda3/envs/scanpyenv"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'docker://raschwaa/scanpyenv-minimal:latest' :
        'raschwaa/scanpyenv-minimal:latest' }"

    input:
    tuple val(study_name), path(study_path)
    val gemma_username
    val gemma_password

    output:
    tuple val(study_name), path("**sample_meta.tsv"), emit: meta

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    python ${projectDir}/bin/get_gemma_meta.py \\
        --study_name ${study_name} \\
        --gemma_username ${gemma_username} \\
        --gemma_password ${gemma_password}
    """

    stub:
    """
    mkdir -p ${study_name}
    echo "sample_id\tsample_name" > ${study_name}/sample_meta.tsv
    """
}
