process LOAD_MASK {
    tag "$study_name"
    // label 'process_single'

    input:
    tuple val(study_name), path(mask_file)
    val use_staging

    output:
    path "message.txt", emit: message

    when:
    task.ext.when == null || task.ext.when

    script:
    def gemma_cmd = use_staging ? "gemma-cli-staging" : "gemma-cli"
    """
    ${gemma_cmd} loadSingleCellData --load-cell-level-characteristics \
        -e ${study_name} \
        -clcFile "${mask_file}" \
        -replaceClc \
        -ignoreSamplesLackingData \
        --data-type NULL \
        -clcName mask \
        -clcDefaultValue false \
        2>> "message.txt"
    """

    stub:
    """
    echo "Loaded Mask for ${study_name}" > message.txt
    """
}
