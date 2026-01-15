process PUBLISH_MULTIQC {
    tag "$study_name"
    // label 'process_single'

    input:
    tuple val(study_name), path(multiqc_html)
    val use_staging
    val version
    val nmads

    output:
    path "message.txt", emit: message

    when:
    task.ext.when == null || task.ext.when

    script:
    def gemma_cmd = use_staging ? "gemma-cli-staging" : "gemma-cli"
    def file_type = use_staging ? "CELL_TYPE_ANNOTATION_PIPELINE_REPORT" : "MULTIQC_REPORT"
    """
    ${gemma_cmd} addMetadataFile \\
        -e ${study_name} \\
        --file-type ${file_type} ${multiqc_html} \\
        --force \\
        --changelog-entry "sc-pipeline-${version} --nmads ${nmads}" \\
        2> message.txt
    """

    stub:
    """
    echo "Published MultiQC for ${study_name}" > message.txt
    """
}
