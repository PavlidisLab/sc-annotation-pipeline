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
    echo "process_samples: ${params.process_samples}" >> params.txt
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
    path "refs/*.h5ad", emit: ref_paths
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
    tag "$query_name"

    publishDir (
        path: "${params.outdir}/${study_name}/predicted_celltypes",
        mode: "copy"
    )

    input:
    tuple val(study_name), val(query_name), val(query_path), val(ref_path)

    output:
    tuple val(study_name), val(query_name), path("${query_name}_*_cell_type.tsv"), emit: celltype_file_channel

    script:
    ref_keys=params.ref_keys.join(" ")
    """
    python $projectDir/bin/scvi_classify.py \\
            --query_path ${query_path} \\
            --ref_path ${ref_path} \\
            --cutoff ${params.cutoff} \\
            --mapping_file ${params.rename_file} \\
            --ref_keys ${ref_keys}

    """

}

process combineCTA {
   publishDir (
        "${params.outdir}/${study_name}/predicted_celltypes", mode: 'copy'
    )
     input:
        tuple val(study_name), val(level),val(query_names), path(combined_celltype_files)

    output:
        tuple val(study_name), path("${study_name}_${level}_combined_celltypes.tsv"), emit: celltype_file_channel

    script:
    """
    # Combine all celltype files into one and only take header from the first file
     # Extract header from first file
    head -n 1 \$(ls ${combined_celltype_files} | head -n 1) > ${study_name}_${level}_combined_celltypes.tsv

    # Append all lines excluding header from all files
    for f in ${combined_celltype_files}; do
        tail -n +2 "\$f" >> ${study_name}_${level}_combined_celltypes.tsv
    done
    """

}

process loadCTA {
    tag "$study_name"

     input:
        tuple val(study_name), path(celltype_file)


    output :
        path "message.txt"



    script:
    def gemma_cmd = params.use_staging ? "gemma-cli-staging" : "gemma-cli"
    def level = celltype_file.getName().split("_")[-3]
    def preferredCtaFlag = (level == "class") ? "-preferredCta" : ""
    """
   ${gemma_cmd} loadSingleCellData -loadCta -e ${study_name} \\
               -ctaFile ${celltype_file} ${preferredCtaFlag} \\
               -ctaName "sc-pipeline-${params.version}-${level}" \\
               -ctaProtocol "sc-pipeline-${params.version}" \\
               -ignoreSamplesLackingData 2> "message.txt"
    """
}

process combineCLC {
   publishDir (
        "${params.outdir}/${study_name}/masks", mode: 'copy'
    )
     input:
        tuple val(study_name), path(combined_mask_files)

    output :
       tuple val(study_name), path("${study_name}_combined_celltype_mask.tsv"), emit: celltype_mask_files

    script:
    """
    # Combine all celltype files into one and only take header from the first file
     # Extract header from first file
    head -n 1 \$(ls ${combined_mask_files} | head -n 1) > ${study_name}_combined_celltype_mask.tsv

    # Append all lines excluding header from all files
    for f in ${combined_mask_files}; do
        tail -n +2 "\$f" >> ${study_name}_combined_celltype_mask.tsv
    done
    """

}

process loadCLC {
    input:
        tuple val(study_name), path(mask_file)

    output:
        path "message.txt"

    script:
    def gemma_cmd = params.use_staging ? "gemma-cli-staging" : "gemma-cli"
    """
    ${gemma_cmd} loadSingleCellData --load-cell-level-characteristics \\
        -e ${study_name} \\
        -clcFile "${mask_file}" \\
        -replaceClc \\
        -ignoreSamplesLackingData \\
        -clcName counts_outlier,genes_outlier,hb_outlier,mito_outlier,predicted_doublet,ribo_outlier,umi_outlier \\
        2>> "message.txt"
    """
}

process getMeta {
    tag "$study_name"
    input:
        tuple val(study_name), path(study_path)

    output:
       tuple val(study_name), path("**sample_meta.tsv"), emit: meta_channel

    script:
    """
    python $projectDir/bin/get_gemma_meta.py --study_name ${study_name} --gemma_username ${params.GEMMA_USERNAME} --gemma_password ${params.GEMMA_PASSWORD}
    """
}


process processQC {
    tag "$query_name"

    publishDir (
        "${params.outdir}/${study_name}/qc", mode: 'copy'
    )

    input:
        tuple val(study_name), val(query_name), path(predicted_meta_files), path(study_path), path(sample_meta)

    output:
    path "**png"
    tuple val(study_name), path("${query_name}/"), emit: qc_channel
    tuple val(study_name), path("${query_name}**mask.tsv"), emit: mask_files


    script:
    """
    python $projectDir/bin/process_QC.py --query_path ${study_path} \\
        --assigned_celltypes_paths ${predicted_meta_files.join(' ')} \\
        --gene_mapping ${params.gene_mapping} \\
        --rename_file ${params.rename_file} \\
        --nmads ${params.nmads} \\
        --sample_meta ${sample_meta} \\
        --organism ${params.organism} \\
        --markers_file ${params.markers_file} \\
        --cell_type_key ${params.cell_type_key}
    """ 
}

process combineQC {
    publishDir (
        "${params.outdir}/${study_name}/combined_qc", mode: 'copy'
    )

    input:
        tuple val(study_name), path(qc_dirs)


    output:
       tuple val(study_name), path("${study_name}/"), emit: qc_dir_combined

    script:
    """
    # combine all qc directories into one directory
    mkdir -p ${study_name}
    for dir in ${qc_dirs}; do
        cp -r \$dir/* "${study_name}/"
    done
    """
}

process runMultiQC {
    tag "$study_name"
    
    publishDir (
        "${params.outdir}/multiqc/${study_name}", mode: 'copy'
    )

    input:
        tuple val(study_name), path(qc_dir)

    output:
        tuple val(study_name), path("**multiqc_report.html"), emit: multiqc_html


    script:
    def use_config_flag = params.process_samples ? "" : "--config new_config.yaml"

    """
    # Combine base config with dynamic title
    cp ${params.multiqc_config} new_config.yaml
    echo 'title: "${study_name} CELLxGENE Census ${params.census_version} cutoff ${params.cutoff} MADs ${params.nmads} "' >> new_config.yaml

    multiqc ${qc_dir} -d ${use_config_flag}
    """
}

process publishMultiQC {

    input:
        tuple val(study_name), path(multiqc_html)

    output:
        path "**message.txt"

    script:
    
    def gemma_cmd = params.use_staging ? "gemma-cli-staging" : "gemma-cli"

    """
    ${gemma_cmd} addMetadataFile \
    -e ${study_name} \
    --file-type ${params.use_staging ? "CELL_TYPE_ANNOTATION_PIPELINE_REPORT" : "MULTIQC_REPORT"} ${multiqc_html} \
    --force \
    --changelog-entry "sc-pipeline-${params.version} --nmads ${params.nmads}" \
    2> message.txt
    """

}

include { DOWNLOAD_STUDIES_SUBWF } from "$projectDir/modules/subworkflows/download_studies.nf"
include { PROCESS_QUERY_SAMPLE } from "$projectDir/modules/processes/process_query_samples.nf"
include { PROCESS_QUERY_COMBINED } from "$projectDir/modules/processes/process_query_combined.nf"

// Workflow definition
workflow {


    DOWNLOAD_STUDIES_SUBWF(params.study_names, params.study_paths)
    DOWNLOAD_STUDIES_SUBWF.out.study_channel.set { study_channel }
    
    // Call the setup process to download the model
    model_path = runSetup(params.organism, params.census_version)

    // If process_samples is true, we will process each query sample separately
    // and use a different process
    def processed_queries
    def raw_queries
    if (params.process_samples) {
        // Split study_channel into individual samples
        expanded_channel = study_channel.flatMap { study_name, study_dir ->
                def results = []
                study_dir.eachDir { dir -> results << [study_name, dir.name, dir.toString()] }
                return results
            }
        // Process each query sample separately
        PROCESS_QUERY_SAMPLE(model_path, expanded_channel)
        raw_queries = PROCESS_QUERY_SAMPLE.out.raw_query
        processed_queries = PROCESS_QUERY_SAMPLE.out.processed_query
    } else {
        // Process each query without subsampling
        PROCESS_QUERY_COMBINED(model_path, study_channel)
        raw_queries = PROCESS_QUERY_COMBINED.out.raw_query
        processed_queries = PROCESS_QUERY_COMBINED.out.processed_query
        
    }
    
    // Get collection names to pull from census
    ref_collections = params.ref_collections.collect { "\"${it}\"" }.join(' ') 

    // Get reference data and save to files
    getCensusAdata(ref_collections)
    getCensusAdata.out.ref_paths.flatten()
    .set { ref_paths }
    
    // Combine the processed queries with the reference paths
    combos_adata = processed_queries.combine(ref_paths)
    // Process each query-reference pair
    rfClassify(combos_adata)

    celltype_files = rfClassify.out.celltype_file_channel
    if (params.process_samples) {
        // Flatten all celltype files for each study/query, then group by study
        celltype_file_channel = celltype_files.flatMap { study_name, query_name, files ->
            files.collect { file -> 
            level = file.getName().split("_")[-3]
            tuple(study_name, level, query_name, file) }
        }.groupTuple(by: [0,1])
         .set{ combined_celltype_files }
        combineCTA(combined_celltype_files)
        predicted_celltypes = combineCTA.out.celltype_file_channel 
    } else {
        // Use the celltype files as they are, flattening the list
        predicted_celltypes = celltype_files.flatMap { study_name, query_name, files ->
            files.collect { file -> tuple(study_name, file) }
        }
    }
    loadCTA(predicted_celltypes)
    predicted_celltypes.groupTuple(by: 0)
    .set { grouped_predicted_celltypes }
 
    getMeta(study_channel)
    meta_channel = getMeta.out.meta_channel


    // need to combine all three of these into one channel with study name, query names, raw query paths, predicted celltype paths, and meta paths
    combined_channel = grouped_predicted_celltypes.combine(raw_queries, by: 0)
    combined_channel = combined_channel.combine(meta_channel, by: 0)

    // put combined channel in the right order if process_samples is true
    if (params.process_samples) {
        combined_channel = combined_channel.map { study_name, predicted_celltypes_paths, query_name, raw_query_path, meta_path ->
            return tuple(study_name, query_name, predicted_celltypes_paths, raw_query_path, meta_path)
        }
    } else {
        combined_channel = combined_channel.map { study_name, predicted_celltypes_paths, query_name, raw_query_path, meta_path ->
            tuple(study_name, query_name, predicted_celltypes_paths, raw_query_path, meta_path)
        }
    }
     

    processQC(combined_channel)

    qc_channel = processQC.out.qc_channel
    mask_files = processQC.out.mask_files
    if (params.process_samples) {
        // If process_samples is true, we will combine the mask files
        // for each study into one file
        // need to combine mask files for each study
        mask_files.groupTuple(by: 0)
        .set{ combined_mask_files }
        combineCLC(combined_mask_files)
        celltype_mask_files = combineCLC.out.celltype_mask_files
    } else {
        // If process_samples is false, we will use the mask files as they are
        celltype_mask_files = mask_files
    } 


    //if (params.mask) {
        // If mask is true, we will load the cell-level characteristics
        //loadCLC(celltype_mask_files)
    //} 


    if (params.process_samples) {
        qc_channel.groupTuple(by: 0)
        .set { qc_dir_channel }
        combineQC(qc_dir_channel)
        multiqc_channel = combineQC.out.qc_dir_combined
    } else {
        // If process_samples is false, we will use the qc_dir as it is
        multiqc_channel = qc_channel
    }


    //// Run MultiQC on the combined qc directory
    if (params.process_samples == false) {
        runMultiQC(multiqc_channel)
        multiqc_channel = runMultiQC.out.multiqc_html
        publishMultiQC(multiqc_channel)
    }
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
