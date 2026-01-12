/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    QC ANALYSIS AND MULTIQC REPORTING
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { GET_META     } from "$projectDir/modules/local/get_meta/main"
include { PROCESS_QC   } from "$projectDir/modules/local/process_qc/main"
include { COMBINE_QC   } from "$projectDir/modules/local/combine_qc/main"
include { COMBINE_CLC  } from "$projectDir/modules/local/combine_clc/main"
include { RUN_MULTIQC  } from "$projectDir/modules/local/run_multiqc/main"
include { PUBLISH_MULTIQC } from "$projectDir/modules/local/publish_multiqc/main"
workflow QC_REPORTING {

    take:
    ch_studies           // channel: [ study_name, study_path ]
    ch_raw               // channel: [ study_name, query_name, raw.h5ad ]
    ch_celltypes         // channel: [ study_name, celltype.tsv ]
    gene_mapping         // path: gene mapping file
    rename_file          // path: cell type renaming file
    nmads                // integer: number of MADs for outlier detection
    organism             // string: organism name
    markers_file         // path: marker genes file
    cell_type_keys       // list: cell type annotation keys
    multiqc_config       // path: MultiQC config file
    census_version       // string: census version
    cutoff               // float: classification cutoff
    process_samples      // boolean: process individual samples
    gemma_username       // string: Gemma username
    gemma_password       // string: Gemma password

    main:
    ch_versions = Channel.empty()

    // Get metadata from Gemma
    GET_META(ch_studies, gemma_username, gemma_password)
    ch_meta = GET_META.out.meta
    ch_versions = ch_versions.mix(GET_META.out.versions.first())

    // Group celltypes by study
    ch_celltypes_grouped = ch_celltypes.groupTuple(by: 0)

    // Combine channels for QC processing
    ch_qc_input = ch_celltypes_grouped
        .combine(ch_raw, by: 0)
        .combine(ch_meta, by: 0)

    // Reorder for PROCESS_QC input
    ch_qc_input = ch_qc_input.map { study_name, celltypes, query_name, raw_path, meta_path ->
        tuple(study_name, query_name, celltypes, raw_path, meta_path)
    }

    // Run QC analysis
    PROCESS_QC(
        ch_qc_input,
        gene_mapping,
        rename_file,
        nmads,
        organism,
        markers_file,
        cell_type_keys
    )
    ch_qc_dirs   = PROCESS_QC.out.qc_dir
    ch_mask_files = PROCESS_QC.out.mask_files
    ch_versions  = ch_versions.mix(PROCESS_QC.out.versions.first())

    if (process_samples) {
        // Combine QC directories
        ch_qc_grouped = ch_qc_dirs.groupTuple(by: 0)
        COMBINE_QC(ch_qc_grouped)
        ch_multiqc_input = COMBINE_QC.out.combined_qc

        // Combine mask files
        ch_mask_grouped = ch_mask_files.groupTuple(by: 0)
        COMBINE_CLC(ch_mask_grouped)
        ch_masks = COMBINE_CLC.out.combined_mask
    }
    else {
        ch_multiqc_input = ch_qc_dirs
        ch_masks = ch_mask_files
    }

    // Run MultiQC (only when not processing samples individually)
    if (!process_samples) {
        RUN_MULTIQC(
            ch_multiqc_input,
            multiqc_config,
            census_version,
            cutoff,
            nmads,
            process_samples
        )

        ch_multiqc = RUN_MULTIQC.out.report
        ch_versions = ch_versions.mix(RUN_MULTIQC.out.versions.first())
    }
    else {
        ch_multiqc = Channel.empty()
    }

    emit:
    qc_dirs   = ch_multiqc_input  // channel: [ study_name, qc_dir ]
    masks     = ch_masks          // channel: [ study_name, mask.tsv ]
    multiqc   = ch_multiqc        // channel: [ study_name, multiqc_report.html ]
    versions  = ch_versions       // channel: [ versions.yml ]
}
