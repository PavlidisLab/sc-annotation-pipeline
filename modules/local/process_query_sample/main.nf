process PROCESS_QUERY_SAMPLE {
    tag "$query_name"
    // label 'process_high'
    errorStrategy { task.exitStatus == 42 ? 'ignore' : 'terminate' }

    conda "/home/rschwartz/anaconda3/envs/scanpyenv"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://raschwaa/census_pipeline:latest' :
        'raschwaa/census_pipeline:latest' }"

    input:
    path model_path
    tuple val(study_name), val(query_name), path(query_path)
    val seed

    output:
    tuple val(study_name), val(query_name), path("${query_name}**.h5ad")    , emit: processed_query
    tuple val(study_name), val(query_name), path("${query_name}_raw.h5ad"), emit: raw_query
    path "versions.yml"                                                   , emit: versions

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

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        scvi-tools: \$(python -c "import scvi; print(scvi.__version__)")
        scanpy: \$(python -c "import scanpy; print(scanpy.__version__)")
    END_VERSIONS
    """

    stub:
    """
    touch ${query_name}.h5ad
    touch ${query_name}_raw.h5ad

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        scvi-tools: \$(python -c "import scvi; print(scvi.__version__)")
        scanpy: \$(python -c "import scanpy; print(scanpy.__version__)")
    END_VERSIONS
    """
}
