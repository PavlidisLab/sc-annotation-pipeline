process COMBINE_QC {
    tag "$study_name"
    // label 'process_single'

    input:
    tuple val(study_name), path(qc_dirs)

    output:
    tuple val(study_name), path("${study_name}/"), emit: combined_qc

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    mkdir -p ${study_name}
    for dir in ${qc_dirs}; do
        cp -r \$dir/* "${study_name}/"
    done
    """

    stub:
    """
    mkdir -p ${study_name}
    touch ${study_name}/qc_combined.tsv
    """
}
