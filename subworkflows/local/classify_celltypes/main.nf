/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CLASSIFY CELL TYPES (RANDOM FOREST OR KNN ON SCVI EMBEDDINGS, OR SEURAT LABEL TRANSFER)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CLASSIFY_SCVI  } from "$projectDir/modules/local/classify_scvi/main"
include { REF_PROCESS_SEURAT   } from "$projectDir/modules/local/ref_process_seurat/main"
include { QUERY_PROCESS_SEURAT } from "$projectDir/modules/local/query_process_seurat/main"
include { PREDICT_SEURAT       } from "$projectDir/modules/local/predict_seurat/main"
include { CLASSIFY_SEURAT      } from "$projectDir/modules/local/classify_seurat/main"
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

    if (params.classifier == 'seurat') {
        // Seurat label transfer: the cutoff applies to the top prediction score, as for rf and knn
        REF_PROCESS_SEURAT(ch_refs)
        QUERY_PROCESS_SEURAT(ch_processed)
        ch_seurat_combos = QUERY_PROCESS_SEURAT.out.query_rds.combine(REF_PROCESS_SEURAT.out.ref_rds)
        PREDICT_SEURAT(ch_seurat_combos, ref_keys)
        CLASSIFY_SEURAT(PREDICT_SEURAT.out.scores, cutoff, mapping_file, ref_keys)
        ch_celltype_files = CLASSIFY_SEURAT.out.celltype_files
    }
    else {
        // Combine processed queries with reference paths
        ch_combos = ch_processed.combine(ch_refs)

        // Classify with the random forest or kNN set by params.classifier
        CLASSIFY_SCVI(ch_combos, cutoff, mapping_file, ref_keys)
        ch_celltype_files = CLASSIFY_SCVI.out.celltype_files
    }

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
