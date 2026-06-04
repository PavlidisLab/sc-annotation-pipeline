process SETUP_SCVI {
    tag "$organism"
    // label 'process_medium'

    conda "/home/rschwartz/anaconda3/envs/scanpyenv"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://raschwaa/scanpyenv-minimal:latest' :
        'raschwaa/scanpyenv-minimal:latest' }"

    input:
    val organism
    val census_version

    output:
    path "scvi-${organism}-${census_version}/", emit: model_path

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    python ${projectDir}/bin/setup.py \\
        --organism ${organism} \\
        --census_version ${census_version} \\
        ${args}
    """

    stub:
    """
    mkdir -p scvi-${organism}-${census_version}
    touch scvi-${organism}-${census_version}/model.pt
    """
}
