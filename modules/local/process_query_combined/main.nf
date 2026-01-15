process PROCESS_QUERY_COMBINED {
    tag "$study_name"
    // label 'process_high'

    conda "/home/rschwartz/anaconda3/envs/scanpyenv"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://raschwaa/census_pipeline:latest' :
        'raschwaa/census_pipeline:latest' }"

    input:
    path model_path
    tuple val(study_name), path(study_path)
    val seed

    output:
    tuple val(study_name), val(study_name), path("${study_name}.h5ad")    , emit: processed_query
    tuple val(study_name), val(study_name), path("${study_name}_raw.h5ad"), emit: raw_query
    path "versions.yml"                                                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    python ${projectDir}/bin/process_query.py \\
        --model_path ${model_path} \\
        --study_path ${study_path} \\
        --study_name ${study_name} \\
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
    touch ${study_name}.h5ad
    touch ${study_name}_raw.h5ad

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        scvi-tools: \$(python -c "import scvi; print(scvi.__version__)")
        scanpy: \$(python -c "import scanpy; print(scanpy.__version__)")
    END_VERSIONS
    """
}
