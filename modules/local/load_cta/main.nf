process LOAD_CTA {
    tag "$study_name GEMMA_CLI_TASK"
    label 'process_single'

    input:
    tuple val(study_name), path(celltype_file)
    val use_staging
    val protocol
    val protocol_description
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
    # Create the protocol if it does not exist yet. The description is read from a quoted heredoc
    # so no character in it is interpreted by the shell. A duplicate name makes addProtocol exit 2
    # with "There is already a protocol object with this name"; only that case is ignored.
    read -r -d '' PROTOCOL_DESCRIPTION <<'PROTOCOL_DESCRIPTION_EOF' || true
${protocol_description}
PROTOCOL_DESCRIPTION_EOF

    set +e
    ${gemma_cmd} addProtocol -name "${protocol}" -description "\$PROTOCOL_DESCRIPTION" > "protocol.txt" 2>&1
    protocol_status=\$?
    set -e
    cat protocol.txt > message.txt
    if [ \$protocol_status -ne 0 ]; then
        if grep -q "There is already a protocol object with this name" protocol.txt; then
            echo "Protocol ${protocol} already exists, continuing" >> message.txt
        else
            echo "addProtocol failed with exit code \$protocol_status" >&2
            tail -n 20 protocol.txt >&2
            exit \$protocol_status
        fi
    fi

    ${gemma_cmd} loadSingleCellData -loadCta -e ${study_name} \\
        -ctaFile ${celltype_file} ${preferred_flag} \\
        -ctaName "${protocol}-${level}" \\
        -ctaProtocol "${protocol}" \\
        -replaceCta \\
        --data-type NULL \\
        -ignoreSamplesLackingData 2>> "message.txt"
    """

    stub:
    """
    echo "Stub: would add protocol ${protocol} if missing" > message.txt
    echo "Loaded CTA for ${study_name}" >> message.txt
    """
}
