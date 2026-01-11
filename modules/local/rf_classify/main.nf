process RF_CLASSIFY {
    tag "$query_name"
    label 'process_medium'

    conda "/home/rschwartz/anaconda3/envs/scanpyenv"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://raschwaa/census_pipeline:latest' :
        'raschwaa/census_pipeline:latest' }"

    input:
    tuple val(study_name), val(query_name), path(query_path), path(ref_path)
    val cutoff
    path mapping_file
    val ref_keys

    output:
    tuple val(study_name), val(query_name), path("${query_name}_*_cell_type.tsv"), emit: celltype_files
    path "versions.yml"                                                          , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def ref_keys_str = ref_keys.join(' ')
    """
    python ${projectDir}/bin/scvi_classify.py \\
        --query_path ${query_path} \\
        --ref_path ${ref_path} \\
        --cutoff ${cutoff} \\
        --mapping_file ${mapping_file} \\
        --ref_keys ${ref_keys_str} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        scikit-learn: \$(python -c "import sklearn; print(sklearn.__version__)")
    END_VERSIONS
    """

    stub:
    """
    touch ${query_name}_subclass_cell_type.tsv
    touch ${query_name}_class_cell_type.tsv
    touch ${query_name}_family_cell_type.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        scikit-learn: \$(python -c "import sklearn; print(sklearn.__version__)")
    END_VERSIONS
    """
}
