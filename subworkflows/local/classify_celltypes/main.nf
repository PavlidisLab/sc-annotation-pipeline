/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CLASSIFY CELL TYPES ON SCVI EMBEDDINGS (RANDOM FOREST OR KNN)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CLASSIFY_SCVI  } from "$projectDir/modules/local/classify_scvi/main"
include { COMBINE_CTA  } from "$projectDir/modules/local/combine_cta/main"

workflow CLASSIFY_CELLTYPES {

    take:
    ch_processed     // channel: [ study_name, query_name, processed.h5ad ]
    ch_refs          // channel: [ ref.h5ad, ... ]
    cutoff           // float: classification probability cutoff
    mapping_file     // path: cell type mapping file
    ref_keys         // list: reference annotation keys
    process_samples  // boolean: whether processing individual samples

    main:

    // Combine processed queries with reference paths
    ch_combos = ch_processed.combine(ch_refs)

    // Classify with the random forest or kNN set by params.classifier
    CLASSIFY_SCVI(ch_combos, cutoff, mapping_file, ref_keys)
    ch_celltype_files = CLASSIFY_SCVI.out.celltype_files

    if (process_samples) {
        // Group and combine cell type files by study and level
        ch_grouped = ch_celltype_files
            .flatMap { study_name, query_name, files ->
                files.collect { file ->
                    def level = file.getName().split("_")[-3]
                    tuple(study_name, level, query_name, file)
                }
            }
            .groupTuple(by: [0, 1])

        COMBINE_CTA(ch_grouped)
        ch_celltypes = COMBINE_CTA.out.combined_celltypes
    }
    else {
        // Flatten files for non-sample processing
        ch_celltypes = ch_celltype_files
            .flatMap { study_name, query_name, files ->
                files.collect { file -> tuple(study_name, file) }
            }
    }

    emit:
    celltypes        = ch_celltypes        // channel: [ study_name, celltype.tsv ]
    celltype_files   = ch_celltype_files   // channel: [ study_name, query_name, [files] ]
}
