/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PROCESS QUERY DATA THROUGH scVI MODEL
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { PROCESS_QUERY_SAMPLE   } from '../../../modules/local/process_query_sample/main'
include { PROCESS_QUERY_COMBINED } from '../../../modules/local/process_query_combined/main'

workflow PROCESS_QUERIES {

    take:
    ch_studies       // channel: [ study_name, study_path ]
    model_path       // path: scVI model directory
    process_samples  // boolean: process individual samples
    seed             // integer: random seed

    main:
    ch_versions = Channel.empty()

    if (process_samples) {
        // Expand channel to individual samples
        ch_expanded = ch_studies.flatMap { study_name, study_dir ->
            def results = []
            study_dir.eachDir { dir ->
                results << [study_name, dir.name, dir.toString()]
            }
            return results
        }

        // Process each sample separately
        PROCESS_QUERY_SAMPLE(model_path, ch_expanded, seed)
        ch_processed = PROCESS_QUERY_SAMPLE.out.processed_query
        ch_raw       = PROCESS_QUERY_SAMPLE.out.raw_query
        ch_versions  = ch_versions.mix(PROCESS_QUERY_SAMPLE.out.versions.first())
    }
    else {
        // Process entire study as one
        PROCESS_QUERY_COMBINED(model_path, ch_studies, seed)
        ch_processed = PROCESS_QUERY_COMBINED.out.processed_query
        ch_raw       = PROCESS_QUERY_COMBINED.out.raw_query
        ch_versions  = ch_versions.mix(PROCESS_QUERY_COMBINED.out.versions.first())
    }

    emit:
    processed = ch_processed  // channel: [ study_name, query_name, processed.h5ad ]
    raw       = ch_raw        // channel: [ study_name, query_name, raw.h5ad ]
    versions  = ch_versions   // channel: [ versions.yml ]
}
