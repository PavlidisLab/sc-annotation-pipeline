process LOAD_CTA {
    tag "$study_name"
    label 'process_single'

    input:
    tuple val(study_name), path(celltype_file)
    val use_staging
    val version
    val preferred_cta_level

    output:
    path "message.txt", emit: message

    when:
    task.ext.when == null || task.ext.when

    script:
    def gemma_cmd = use_staging ? "gemma-cli-staging" : "gemma-cli"
    def level = celltype_file.getName().split("_")[-3]
    def preferred_flag = (level == preferred_cta_level) ? "-preferredCta" : ""
    """
    ${gemma_cmd} loadSingleCellData -loadCta -e ${study_name} \\
        -ctaFile ${celltype_file} ${preferred_flag} \\
        -ctaName "sc-pipeline-${version}-${level}" \\
        -ctaProtocol "sc-pipeline-${version}" \\
        -replaceCta \\
        --data-type NULL \\
        -ignoreSamplesLackingData 2> "message.txt"
    """

    stub:
    """
    echo "Loaded CTA for ${study_name}" > message.txt
    """
}
