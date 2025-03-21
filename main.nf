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
    echo "studies_dir: ${params.studies_dir}" >> params.txt
    echo "subsample ref: ${params.subsample_ref}" >> params.txt
    echo "ref collections: ${params.ref_collections}" >> params.txt
    """
}



 process downloadStudies {
   // publishDir "${params.outdir}/studies", mode: 'copy'

    input:
        val study_name

    output:
        tuple val(study_name), path("${study_name}/"), emit: study_channel


     script:

     """
     gemma-cli-sc getSingleCellDataMatrix -e $study_name --format mex --scale-type count --use-ensembl-ids -o $study_name
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
    path "${study_name}.h5ad", emit: processed_query

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
    val organism
    val census_version
    val subsample_ref
    val ref_collections
    val organ

    output:
    path "refs/*.h5ad", emit: ref_paths_adata
    path "refs/ref_cell_info.tsv"

    script:
    """
    # Run the python script to generate the files
    python $projectDir/bin/get_census_adata.py \\
        --organism ${organism} \\
        --organ ${organ} \\
        --census_version ${census_version} \\
        --subsample_ref ${subsample_ref} \\
        --ref_collections ${ref_collections} \\
        --rename_file ${params.rename_file} \\
        --seed ${params.seed} \\ 
        --restricted_celltypes ${params.restricted_celltypes}

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
    tuple val(query_path), val(ref_path)

    output:
    path "${query_path.getName().toString().replace(".h5ad","")}/${query_path.getName().toString().replace(".h5ad","")}_predicted_celltype.tsv"

    script:
    """
    python $projectDir/bin/scvi_classify.py --query_path ${query_path} --ref_path ${ref_path} --cutoff ${params.cutoff}    
 
    """

}

// process loadResults {
    // input:
        // tuple path(celltype_file), val(study_name), val(target_platform)


    //output :
        // path message.txt

    // script:
    // """
    // gemma-cli-sc loadSingleCellData -e ${study_name} -a ${params.target_platform} --preferred-quantitation-type
    // """

// Workflow definition
workflow {


    // Get query names from file (including region)
    study_names = Channel.fromPath(params.study_names).flatMap { file ->
        // Read the file, split by lines, and trim any extra spaces
        file.readLines().collect { it.trim() }
    }

    downloadStudies(study_names)
    downloadStudies.out.study_channel.set { study_channel }
    study_channel.view()

    study_channel.view()
    // Call the setup process to download the model
    model_path = runSetup(params.organism, params.census_version)

    // Process each query by relabeling, subsampling, and passing through scvi model
    processed_queries_adata = processQuery(model_path, study_channel) 
     
    // Get collection names to pull from census
    ref_collections = params.ref_collections.collect { "\"${it}\"" }.join(' ')
    
    // Get reference data and save to files
    getCensusAdata(params.organism, params.census_version, params.subsample_ref, ref_collections, params.organ)
    getCensusAdata.out.ref_paths_adata.flatten()
    .set { ref_paths_adata }
    
    // Combine the processed queries with the reference paths
    combos_adata = processed_queries_adata.combine(ref_paths_adata)
    
    // Process each query-reference pair
    rfClassify(combos_adata)

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
