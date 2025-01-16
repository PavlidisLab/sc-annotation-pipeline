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
// process parseJsonMeta {
//     input:
//         path study_json_file
//     output:
//         path study_meta_file
//     script:
//         """
//         python $projectDir/bin/parse_json.py --json_file ${study_json_file}
//         """
// }
//


// process parseStudies {
    // input:
        // path study_meta_file
    // output:
        // tuple val(study_name), val(mapped_organism)

    // script:


    // """
    // study_name = study_meta_file.split("_")[0]
    // readLines("${study_meta_file}").first() { line ->
    //
        // def organism = line.split(" ")[1]
        //
        //
    //
    //  
        // # Output the tuple
        // echo "$study_name $organism
    // done
    // """
// }


// process getStudies {

    // input:
        // val study_name, val organism

    // output:
        // 
        // path(${params.studies_dir}/${experiment}")
        //

    // script:

    // """
    // gemma-cli-sc getSingleCellDataMatrix -e ${study_name} \\
    // --format mex --scale-type count --use-ensembl-ids \\
    // -o /space/scratch/gemma-single-cell-data-ensembl-id/${organism}/${study_name}
    // """
// }

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
    tuple val(study_path), val(study_name)

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

    publishDir (
       path: "${params.outdir}",
        mode: "copy"
    )

    input:
    val organism
    val census_version
    val subsample_ref
    val ref_collections

    output:
    path "refs/*.h5ad", emit: ref_paths_adata

    script:
    """
    # Run the python script to generate the files
    python $projectDir/bin/get_census_adata.py \\
        --organism ${organism} \\
        --census_version ${census_version} \\
        --subsample_ref ${subsample_ref} \\
        --ref_collections ${ref_collections} \\
        --seed ${params.seed}

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
        // path "*.tsv"
        // val experiment
        // val target_platform

    //output :
        // path message.txt

    // script:
    // """
    // gemma-cli loadSingleCellData -e <experiment ID> -p <target platform>
    // """

// Workflow definition
workflow {

    Channel
        .fromPath("${params.studies_dir}/*", type: 'dir')
        .set { study_paths }

    // Get query names from file (including region)
    study_paths = study_paths.map{ study_path -> 
        def study_name = study_path.toString().split('/')[-1]
        [study_path, study_name]
    }

    // study_meta = parseStudies(params.study_meta_file)

    // study_channel = getStudies(study_meta)

    // combined_study_channel = study_channel.map{ study_path -> 
        // def study_name = study_path.toString().split('/')[-1]
        // def organism = study_path.toString().split('/')[-2]
        // [study_path, study_name, organism]
    // }


    // Call the setup process to download the model
    model_path = runSetup(params.organism, params.census_version)

    // Process each query by relabeling, subsampling, and passing through scvi model
    processed_queries_adata = processQuery(model_path, study_paths) 
     
    // Get collection names to pull from census
    ref_collections = params.ref_collections.collect { "\"${it}\"" }.join(' ')
    
    // Get reference data and save to files
    getCensusAdata(params.organism, params.census_version, params.subsample_ref, ref_collections)
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
