process LOAD_CLC {
    tag "$study_name GEMMA_CLI_TASK"
    label 'process_single'

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
    # Gemma names the characteristics in the order their category first appears in the file, so the names are
    # read from the file rather than listed here; a fixed list scrambled the names whenever the two orders differed.
    names=\$(awk -F'\\t' 'NR == 1 { for (i = 1; i <= NF; i++) if (\$i == "category") c = i; next } !seen[\$c]++ { printf "%s%s", (n++ ? "," : ""), \$c }' "${clc_file}")
    defaults=\$(echo "\$names" | sed 's/[^,]*/false/g')
    ${gemma_cmd} loadSingleCellData --load-cell-level-characteristics \\
        -e ${study_name} \\
        -clcFile "${clc_file}" \\
        -replaceClc \\
        -ignoreSamplesLackingData \\
        --data-type NULL \\
        -clcName "\$names" \\
        -clcDefaultValue "\$defaults" \\
        2>> "message.txt"
    """

    stub:
    """
    echo "Loaded CLC for ${study_name}" > message.txt
    """
}
