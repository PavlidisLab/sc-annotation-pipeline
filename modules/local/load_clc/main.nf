process LOAD_CLC {
    tag "$study_name"
    // label 'process_single'

    input:
    tuple val(study_name), path(clc_file)
    val use_staging

    output:
    path "message.txt", emit: message

    when:
    task.ext.when == null || task.ext.when

    script:
    def gemma_cmd = use_staging ? "gemma-cli-staging" : "gemma-cli"
    """
    ${gemma_cmd} loadSingleCellData --load-cell-level-characteristics \\
        -e ${study_name} \\
        -clcFile "${clc_file}" \\
        -replaceClc \\
        -ignoreSamplesLackingData \\
        --data-type NULL \\
        -clcName counts_outlier,genes_outlier,mito_outlier,predicted_doublet,non_outlier,umi_outlier \\
        -clcDefaultValue false,false,false,false,false,false \\
        2>> "message.txt"
    """

    stub:
    """
    echo "Loaded CLC for ${study_name}" > message.txt
    """
}
