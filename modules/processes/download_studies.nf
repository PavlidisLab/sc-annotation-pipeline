process DOWNLOAD_STUDIES {
    publishDir "${params.outdir}/mex", mode: 'copy'

    tag "$study_name"

    input:
        val study_name

    output:
        tuple val(study_name), path("${study_name}/")

    script:
    """
    gemma-cli-staging getSingleCellDataMatrix -e $study_name --format mex --scale-type count --use-ensembl-ids -o $study_name
    """
}