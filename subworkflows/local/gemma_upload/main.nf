/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    UPLOAD RESULTS TO GEMMA DATABASE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { LOAD_CTA        } from "$projectDir/modules/local/load_cta/main"
include { LOAD_CLC        } from "$projectDir/modules/local/load_clc/main"
include { PUBLISH_MULTIQC } from "$projectDir/modules/local/publish_multiqc/main"

workflow GEMMA_UPLOAD {

    take:
    ch_celltypes          // channel: [ study_name, celltype.tsv ]
    ch_masks              // channel: [ study_name, mask.tsv ]
    ch_multiqc            // channel: [ study_name, multiqc_report.html ]
    use_staging           // boolean: use Gemma staging environment
    version               // string: pipeline version
    preferred_cta_level   // string: preferred CTA level
    nmads                 // integer: MADs used
    upload_cta            // boolean: upload cell type annotations
    upload_clc            // boolean: upload cell-level characteristics
    upload_multiqc        // boolean: upload MultiQC report

    main:
    ch_messages = Channel.empty()

    // Upload cell type annotations
    if (upload_cta) {
        LOAD_CTA(ch_celltypes, use_staging, version, preferred_cta_level)
        ch_messages = ch_messages.mix(LOAD_CTA.out.message)
    }

    // Upload cell-level characteristics (outlier masks)
    if (upload_clc) {
        LOAD_CLC(ch_masks, use_staging)
        ch_messages = ch_messages.mix(LOAD_CLC.out.message)
    }

    // Upload MultiQC reports
    if (upload_multiqc) {
        PUBLISH_MULTIQC(ch_multiqc, use_staging, version, nmads)
        ch_messages = ch_messages.mix(PUBLISH_MULTIQC.out.message)
    }

    emit:
    messages = ch_messages  // channel: [ message.txt ]
}
