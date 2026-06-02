/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { INPUT_CHECK        } from "$projectDir/subworkflows/local/input_check/main"
include { PREPARE_REFERENCE  } from "$projectDir/subworkflows/local/prepare_reference/main"
include { PROCESS_QUERIES    } from "$projectDir/subworkflows/local/process_queries/main"
include { CLASSIFY_CELLTYPES } from "$projectDir/subworkflows/local/classify_celltypes/main"
include { QC_REPORTING       } from "$projectDir/subworkflows/local/qc_reporting/main"
include { GEMMA_UPLOAD       } from "$projectDir/subworkflows/local/gemma_upload/main"

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow SCANNOTATE {

    take:
    input         // path: samplesheet CSV (optional)
    study_names   // string: space-separated or file (optional, legacy)
    study_paths   // string: space-separated or file (optional, legacy)

    main:
    ch_versions = Channel.empty()

    //
    // SUBWORKFLOW: Validate inputs and download/prepare studies
    //
    INPUT_CHECK(
        input,
        study_names,
        study_paths,
        params.use_staging
    )
    ch_studies = INPUT_CHECK.out.studies
    ch_versions = ch_versions.mix(INPUT_CHECK.out.versions)

    //
    // SUBWORKFLOW: Prepare reference data (scVI model + Census data)
    //
    PREPARE_REFERENCE(
        params.organism,
        params.census_version,
        params.ref_collections,
        params.organ,
        params.subsample_ref,
        file(params.rename_file),
        params.seed,
        params.original_celltype_columns ? file(params.original_celltype_columns) : [],
        params.author_annotations_path ? file(params.author_annotations_path) : []
    )
    ch_model = PREPARE_REFERENCE.out.model_path
    ch_refs  = PREPARE_REFERENCE.out.ref_paths
    ch_versions = ch_versions.mix(PREPARE_REFERENCE.out.versions)

    //
    // SUBWORKFLOW: Process query data through scVI model
    //
    PROCESS_QUERIES(
        ch_studies,
        ch_model,
        params.process_samples,
        params.seed
    )
    ch_processed = PROCESS_QUERIES.out.processed
    ch_raw       = PROCESS_QUERIES.out.raw
    ch_versions  = ch_versions.mix(PROCESS_QUERIES.out.versions)

    //
    // SUBWORKFLOW: Classify cell types using random forest
    //
    CLASSIFY_CELLTYPES(
        ch_processed,
        ch_refs,
        params.cutoff,
        file(params.rename_file),
        params.ref_keys,
        params.process_samples
    )
    ch_celltypes = CLASSIFY_CELLTYPES.out.celltypes
    ch_versions  = ch_versions.mix(CLASSIFY_CELLTYPES.out.versions)

    //
    // SUBWORKFLOW: QC analysis and MultiQC reporting
    //
    QC_REPORTING(
        ch_studies,
        ch_raw,
        ch_celltypes,
        file(params.gene_mapping),
        file(params.rename_file),
        params.nmads,
        params.organism,
        file(params.markers_file),
        params.ref_keys,
        file(params.multiqc_config),
        params.census_version,
        params.cutoff,
        params.process_samples,
        params.GEMMA_USERNAME,
        params.GEMMA_PASSWORD
    )
    ch_clc     = QC_REPORTING.out.clc
    ch_masks   = QC_REPORTING.out.masks
    ch_multiqc = QC_REPORTING.out.multiqc
    ch_versions = ch_versions.mix(QC_REPORTING.out.versions)

    //
    // SUBWORKFLOW: Upload results to Gemma (optional)
    //
    // Master switch: when --upload_gemma is false the whole upload subworkflow is
    // skipped, so a run needs no Gemma write access. The granular upload_cta/clc/mask/
    // multiqc flags still select which artifacts to upload when on. (Note: fetching study
    // data and metadata from Gemma still requires GEMMA_USERNAME/GEMMA_PASSWORD.)
    //
    ch_messages = Channel.empty()
    if (params.upload_gemma) {
        GEMMA_UPLOAD(
            ch_celltypes,
            ch_clc,
            ch_masks,
            ch_multiqc,
            params.use_staging,
            params.version,
            params.preferredCtaLevel,
            params.nmads,
            params.upload_cta ?: false,
            params.upload_clc ?: false,
            params.upload_mask ?: false,
            params.upload_multiqc ?: false,
            params.process_samples
        )
        ch_messages = GEMMA_UPLOAD.out.messages
    }


    emit:
    celltypes = ch_celltypes
    masks     = ch_masks
    multiqc   = ch_multiqc
    versions  = ch_versions
    messages  = ch_messages
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
