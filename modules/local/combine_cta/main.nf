process COMBINE_CTA {
    tag "$study_name GEMMA_CLI_TASK"
    // label 'process_single'

    input:
    tuple val(study_name), val(level), val(query_names), path(celltype_files)

    output:
    tuple val(study_name), path("${study_name}_${level}_combined_celltypes.tsv"), emit: combined_celltypes

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Extract header from first file
    head -n 1 \$(ls ${celltype_files} | head -n 1) > ${study_name}_${level}_combined_celltypes.tsv

    # Append all lines excluding header from all files
    for f in ${celltype_files}; do
        tail -n +2 "\$f" >> ${study_name}_${level}_combined_celltypes.tsv
    done
    """

    stub:
    """
    echo "cell_id\tcell_type\tprobability" > ${study_name}_${level}_combined_celltypes.tsv
    """
}
