process SETUP_SCVI {
    tag "$organism"
    label 'process_medium'

    conda "/home/rschwartz/anaconda3/envs/scanpyenv"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'docker://raschwaa/census-minimal:latest' :
        'raschwaa/census-minimal:latest' }"

    // model.pt is fully keyed by organism + census_version, so cache it across work-dirs/resumes
    storeDir "${projectDir}/.cache/scvi_model"

    input:
    val organism
    val census_version

    output:
    path "scvi-${organism}-${census_version}/model.pt", emit: model_file

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
