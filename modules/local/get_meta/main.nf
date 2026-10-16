process GET_META {
    tag "$study_name"
    label 'process_single'

    conda "/home/rschwartz/anaconda3/envs/scanpyenv"
    container "${ workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container ?
        'docker://raschwaa/scanpyenv-minimal:latest' :
        'raschwaa/scanpyenv-minimal:latest' }"

    input:
    tuple val(study_name), path(study_path)
    val use_gemma

    output:
    tuple val(study_name), path("**sample_meta.tsv"), emit: meta

    when:
    task.ext.when == null || task.ext.when

    script:
    if (use_gemma)
        // GEMMA_USERNAME / GEMMA_PASSWORD are read by the script from the task environment
        // (inherited from the submitting shell), never passed as arguments.
        """
        python ${projectDir}/bin/get_gemma_meta.py \\
            --study_name ${study_name} \\
            ${params.use_staging ? '--use_staging' : ''}
        """
    else
        // No matching Gemma study to fetch metadata for (--use_gemma false); emit an
        // empty stub so the downstream merge in qc_utils.py's read_query() is a no-op
        // left join instead of failing outright. Needs >=2 tab-separated columns or
        // pandas' `sep=None` delimiter sniffer misparses a single-column header.
        """
        echo -e "sample_id\tsample_name" > ${study_name}_sample_meta.tsv
        """

    stub:
    """
    mkdir -p ${study_name}
    echo "sample_id\tsample_name" > ${study_name}/sample_meta.tsv
    """
}
