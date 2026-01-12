process PROCESS_QC {
    tag "$query_name"
    label 'process_medium'

    conda "/home/rschwartz/anaconda3/envs/scanpyenv"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'docker://raschwaa/census_pipeline:latest' :
        'raschwaa/census_pipeline:latest' }"

    input:
    tuple val(study_name), val(query_name), path(predicted_meta_files), path(study_path), path(sample_meta)
    path gene_mapping
    path rename_file
    val nmads
    val organism
    path markers_file
    val cell_type_keys

    output:
    path "**png"                                    , emit: plots
    tuple val(study_name), path("${query_name}/")   , emit: qc_dir
    tuple val(study_name), path("${query_name}*clc.tsv"), emit: mask_files
    path "versions.yml"                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def cell_type_keys_str = cell_type_keys.join(' ')
    def predicted_files_str = predicted_meta_files.join(' ')
    """
    python ${projectDir}/bin/process_QC.py \\
        --query_path ${study_path} \\
        --assigned_celltypes_paths ${predicted_files_str} \\
        --gene_mapping ${gene_mapping} \\
        --rename_file ${rename_file} \\
        --nmads ${nmads} \\
        --sample_meta ${sample_meta} \\
        --organism ${organism} \\
        --markers_file ${markers_file} \\
        --cell_type_keys ${cell_type_keys_str} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        scanpy: \$(python -c "import scanpy; print(scanpy.__version__)")
    END_VERSIONS
    """

    stub:
    """
    mkdir -p ${query_name}
    touch ${query_name}/qc_metrics.tsv
    touch ${query_name}_mask.tsv
    touch ${query_name}_plot.png

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        scanpy: \$(python -c "import scanpy; print(scanpy.__version__)")
    END_VERSIONS
    """
}
