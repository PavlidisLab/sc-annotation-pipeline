process SETUP_SCVI {
    tag "$organism"
    label 'process_medium'

    conda "/home/rschwartz/anaconda3/envs/scanpyenv"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://raschwaa/census_pipeline:latest' :
        'raschwaa/census_pipeline:latest' }"

    input:
    val organism
    val census_version

    output:
    path "scvi-${organism}-${census_version}/", emit: model_path
    path "versions.yml"                       , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    python ${projectDir}/bin/setup.py \\
        --organism ${organism} \\
        --census_version ${census_version} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        scvi-tools: \$(python -c "import scvi; print(scvi.__version__)")
    END_VERSIONS
    """

    stub:
    """
    mkdir -p scvi-${organism}-${census_version}
    touch scvi-${organism}-${census_version}/model.pt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        scvi-tools: \$(python -c "import scvi; print(scvi.__version__)")
    END_VERSIONS
    """
}
