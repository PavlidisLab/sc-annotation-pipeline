process COMBINE_MASK {
    tag "$study_name"
    // label 'process_single'

    input:
    tuple val(study_name), path(mask_files)

    output:
    tuple val(study_name), path("${study_name}_combined_mask.tsv"), emit: combined_mask

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    # Extract header from first file
    head -n 1 \$(ls ${mask_files} | head -n 1) > ${study_name}_combined_mask.tsv

    # Append all lines excluding header from all files
    for f in ${mask_files}; do
        tail -n +2 "\$f" >> ${study_name}_combined_mask.tsv
    done
    """

    stub:
    """
    echo -e "sample_id\tcell_id\tcategory\tvalue" > ${study_name}_combined_mask.tsv
    """
}
