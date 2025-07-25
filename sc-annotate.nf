#!/usr/bin/env nextflow
process save_params_to_file {

    publishDir (
        "${params.outdir}",
        mode: "copy"
    )

    output:
    file "params.txt"

    script:


    """
    echo "organism: ${params.organism}" > params.txt
    echo "census_version: ${params.census_version}" >> params.txt 
    echo "outdir: ${params.outdir}" >> params.txt
    echo "study_names: ${params.study_names}" >> params.txt
    echo "subsample ref: ${params.subsample_ref}" >> params.txt
    echo "ref collections: ${params.ref_collections}" >> params.txt
    echo "rename file: ${params.rename_file}" >> params.txt
    echo "original celltype columns: ${params.original_celltype_columns}" >> params.txt
    echo "author annotations path: ${params.author_annotations_path}" >> params.txt
    echo "seed: ${params.seed}" >> params.txt
    echo "gene mapping: ${params.gene_mapping}" >> params.txt
    echo "markers file: ${params.markers_file}" >> params.txt
    echo "nmads: ${params.nmads}" >> params.txt
    echo "cutoff: ${params.cutoff}" >> params.txt
    echo "multiqc config: ${params.multiqc_config}" >> params.txt
    echo "version: ${params.version}" >> params.txt

    """
}

process runSetup {
    //conda '/home/rschwartz/anaconda3/envs/scanpyenv'

    input:
    val organism
    val census_version

    output:
    path "scvi-${params.organism}-${census_version}/"

    script:
    """
    python $projectDir/bin/setup.py --organism ${organism} --census_version ${census_version}
    """
}

process processQuery {
  //  conda '/home/rschwartz/anaconda3/envs/scanpyenv'

    input:
    val model_path
    tuple val(study_name), path(study_path)

    output:
    tuple val("${study_name}"), path("${study_name}.h5ad"), emit: processed_query
    tuple val("${study_name}"), path("${study_name}_raw.h5ad"), emit: raw_query

        
    script:


    """

    python $projectDir/bin/process_query.py \\
                            --model_path ${model_path} \\
                            --study_path ${study_path} \\
                            --study_name ${study_name} \\
                            --seed ${params.seed}
    """

}

process getCensusAdata {
  //  conda '/home/rschwartz/anaconda3/envs/scanpyenv'

    input:
    val ref_collections

    output:
    path "refs/*.h5ad", emit: ref_paths_adata
    //path "**ref_cell_info.tsv"

    script:
    """
    # Run the python script to generate the files
    python $projectDir/bin/get_census_adata.py \\
        --organism ${params.organism} \\
        --organ ${params.organ} \\
        --census_version ${params.census_version} \\
        --subsample ${params.subsample_ref} \\
        --ref_collections ${ref_collections} \\
        --rename_file ${params.rename_file} \\
        --seed ${params.seed} \\
        ${params.original_celltype_columns ? "--original_celltype_columns ${params.original_celltype_columns}" : ''} \\
        ${params.author_annotations_path ? "--author_annotations_path ${params.author_annotations_path}" : ''}
        

    # After running the python script, all .h5ad files will be saved in the refs/ directory inside a work directory
    """
}


process rfClassify{
  //  conda '/home/rschwartz/anaconda3/envs/scanpyenv'

    publishDir (
        path: "${params.outdir}",
        mode: "copy"
    )

    input:
    tuple val(study_name), val(query_path), val(ref_path)

    output:
    tuple val{study_name}, path("${study_name}/${study_name}_predicted_celltype.tsv"), emit : celltype_file_channel

    script:
    """
    python $projectDir/bin/scvi_classify.py --query_path ${query_path} --ref_path ${ref_path} --cutoff ${params.cutoff}    
 
    """

}

process loadResults {
    publishDir (
        "${params.outdir}/${study_name}", mode: 'copy'
    )
     input:
        tuple val(study_name), path(celltype_file)


    output :
        path "message.txt"


    """

    gemma-cli loadSingleCellData -loadCta -e ${study_name} \\
               -ctaFile ${celltype_file} -preferredCta \\
               -ctaName "sc-pipeline-${params.version}" \\
               -ctaProtocol "sc-pipeline-${params.version}" 2> "message.txt" 
    """
}

process getMeta {
    //publishDir (
        //"${params.outdir}/${study_name}", mode: 'copy'
    //)

    input:
        tuple val(study_name), path(study_path)

    output:
       tuple val(study_name), path("**sample_meta.tsv"), emit: meta_channel

    script:
    """
    python $projectDir/bin/get_gemma_meta.py --study_name ${study_name} --gemma_username ${params.GEMMA_USERNAME} --gemma_password ${params.GEMMA_PASSWORD}
    """
}


process plotQC {
   // conda '/home/rschwartz/anaconda3/envs/scanpyenv'
    
    //publishDir (
        //"${params.outdir}/${study_name}/qc_results", mode: 'copy'
    //)

    input:
        tuple val(study_name), path(predicted_meta), path(study_path), path(sample_meta)

    output:
    path "**png"
    tuple val(study_name), path("${study_name}/"), emit: qc_channel


    script:
    """
    python $projectDir/bin/plot_QC.py --query_path ${study_path} \\
        --assigned_celltypes_path ${predicted_meta} \\
        --gene_mapping ${params.gene_mapping} \\
        --rename_file ${params.rename_file} \\
        --nmads ${params.nmads} \\
        --sample_meta ${sample_meta} \\
        --organism ${params.organism} \\
        --markers_file ${params.markers_file}
    """ 
}


process runMultiQC {
    publishDir (
        "${params.outdir}/multiqc/${study_name}", mode: 'copy'
    )

    input:
        tuple val(study_name), path(qc_dir)

    output:
        tuple val(study_name), path("multiqc_report.html"), emit: multiqc_html

    script:
    """
    multiqc ${qc_dir} -d --config ${params.multiqc_config}
    """
}

process publishMultiQC {
    publishDir (
        "${params.outdir}/multiqc/${study_name}", mode: 'copy'
    )

    input:
        tuple val(study_name), path(multiqc_html)

    output:
        path "**message.txt"

    script:
    """
    gemma-cli addMetadataFile -e ${study_name} --file-type MULTIQC_REPORT ${multiqc_html} --force --changelog-entry "sc-pipeline-${params.version} --nmads ${params.nmads}" 2> "message.txt"
    """
}

include { DOWNLOAD_STUDIES_SUBWF } from "${projectDir}/modules/subworkflows/download_studies.nf"

// Workflow definition
workflow {


    DOWNLOAD_STUDIES_SUBWF(params.study_names, params.studies_path)

    DOWNLOAD_STUDIES_SUBWF.out.study_channel.set { study_channel }
    
    // Call the setup process to download the model
    model_path = runSetup(params.organism, params.census_version)

    // Process each query by relabeling, subsampling, and passing through scvi model
    processQuery(model_path, study_channel) 
    processed_queries_adata = processQuery.out.processed_query 
    // Get collection names to pull from census
    ref_collections = params.ref_collections.collect { "\"${it}\"" }.join(' ') 

    // Get reference data and save to files
    getCensusAdata(ref_collections)
    getCensusAdata.out.ref_paths_adata.flatten()
    .set { ref_paths_adata }
    
    // Combine the processed queries with the reference paths
    combos_adata = processed_queries_adata.combine(ref_paths_adata)
    // Process each query-reference pair
    rfClassify(combos_adata)

    celltype_files = rfClassify.out.celltype_file_channel

    raw_queries = processQuery.out.raw_query
    celltype_files.join(raw_queries, by: 0)
    .set{qc_channel}

    getMeta(study_channel)
    meta_channel = getMeta.out.meta_channel

    qc_channel.join(meta_channel, by: 0)
    .set { qc_channel_with_meta }
    plotQC(qc_channel_with_meta)
    multiqc_channel = plotQC.out.qc_channel

    runMultiQC(multiqc_channel)

    loadResults(celltype_files)

    multiqc_channel = runMultiQC.out.multiqc_html

    publishMultiQC(multiqc_channel)

    save_params_to_file()
}

workflow.onComplete {
    println "Successfully completed"
    println ( workflow.success ? 
    """
    ===============================================================================
    Pipeline execution summary
    -------------------------------------------------------------------------------

    Run as      : ${workflow.commandLine}
    Started at  : ${workflow.start}
    Completed at: ${workflow.complete}
    Duration    : ${workflow.duration}
    Success     : ${workflow.success}
    workDir     : ${workflow.workDir}
    Config files: ${workflow.configFiles}
    exit status : ${workflow.exitStatus}
    version : ${params.version}
    outdir : ${params.outdir}

    --------------------------------------------------------------------------------
    ================================================================================
    """.stripIndent() : """
    Failed: ${workflow.errorReport}
    exit status : ${workflow.exitStatus}
    """.stripIndent()
    )
}

workflow.onError = {
println "Error: something went wrong, check the pipeline log at '.nextflow.log"
}
