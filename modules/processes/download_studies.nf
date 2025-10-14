process DOWNLOAD_STUDIES {
    //publishDir "${params.cacheDir}/mex", mode: 'copy'

    tag "$study_name"

    input:
        val study_name

    output:
        tuple val(study_name), path("${study_name}/")

    script:
    def gemma_cmd = params.use_staging ? "gemma-cli-staging" : "gemma-cli"
    """
    ${gemma_cmd} getSingleCellDataMatrix -e $study_name --format mex --scale-type count --use-ensembl-ids -o $study_name
    """
}