process PROCESS_QUERY_SAMPLE {
    tag "$query_name"
    label 'process_high'
    errorStrategy { task.exitStatus == 42 ? 'ignore' : 'terminate' }

    conda "/home/rschwartz/anaconda3/envs/scanpyenv"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'docker://raschwaa/scanpyenv-minimal:latest' :
        'raschwaa/scanpyenv-minimal:latest' }"

    input:
    path model_path
    tuple val(study_name), val(query_name), path(query_path)
    val seed

    output:
    tuple val(study_name), val(query_name), path("${query_name}**.h5ad")    , emit: processed_query
    tuple val(study_name), val(query_name), path("${query_name}_raw.h5ad"), emit: raw_query

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    python ${projectDir}/bin/process_query_samples.py \\
        --model_path ${model_path} \\
        --query_name ${query_name} \\
        --query_path ${query_path} \\
        --seed ${seed} \\
        ${args}
    """

    stub:
    """
    touch ${query_name}.h5ad
    touch ${query_name}_raw.h5ad
    """
}
