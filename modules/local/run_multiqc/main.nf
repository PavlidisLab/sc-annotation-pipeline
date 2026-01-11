process RUN_MULTIQC {
    tag "$study_name"
    label 'process_single'

    conda "/home/rschwartz/anaconda3/envs/scanpyenv"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://raschwaa/census_pipeline:latest' :
        'raschwaa/census_pipeline:latest' }"

    input:
    tuple val(study_name), path(qc_dir)
    path multiqc_config
    val census_version
    val cutoff
    val nmads
    val process_samples

    output:
    tuple val(study_name), path("**multiqc_report.html"), emit: report
    path "versions.yml"                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def use_config_flag = process_samples ? "" : "--config new_config.yaml"
    """
    # Combine base config with dynamic title
    cp ${multiqc_config} new_config.yaml
    echo 'title: "${study_name} CELLxGENE Census ${census_version} cutoff ${cutoff} MADs ${nmads}"' >> new_config.yaml

    multiqc ${qc_dir} -d ${use_config_flag} ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$(multiqc --version | sed 's/multiqc, version //g')
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${study_name}
    touch ${study_name}/multiqc_report.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        multiqc: \$(multiqc --version | sed 's/multiqc, version //g')
    END_VERSIONS
    """
}
