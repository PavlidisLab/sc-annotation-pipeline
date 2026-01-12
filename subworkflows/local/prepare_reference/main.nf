/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PREPARE REFERENCE DATA (scVI model and Census AnnData)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { SETUP_SCVI       } from "$projectDir/modules/local/setup_scvi/main"
include { GET_CENSUS_ADATA } from "$projectDir/modules/local/get_census_adata/main"

workflow PREPARE_REFERENCE {

    take:
    organism                  // string: organism name
    census_version           // string: CellxGene census version
    ref_collections          // list: reference collection names
    organ                    // string: organ name
    subsample_ref            // integer: cells per cell type to subsample
    rename_file              // path: cell type renaming file
    seed                     // integer: random seed
    original_celltype_columns // path: optional original celltype columns file
    author_annotations_path   // path: optional author annotations directory

    main:
    ch_versions = Channel.empty()

    // Download scVI model
    SETUP_SCVI(organism, census_version)
    ch_model = SETUP_SCVI.out.model_path
    ch_versions = ch_versions.mix(SETUP_SCVI.out.versions)

    // Get reference data from CellxGene Census
    ref_collections_str = ref_collections.collect { "\"${it}\"" }.join(' ')

    GET_CENSUS_ADATA(
        ref_collections_str,
        organism,
        organ,
        census_version,
        subsample_ref,
        rename_file,
        seed,
        original_celltype_columns ?: [],
        author_annotations_path ?: []
    )
    ch_refs = GET_CENSUS_ADATA.out.ref_paths.flatten()
    ch_versions = ch_versions.mix(GET_CENSUS_ADATA.out.versions)

    emit:
    model_path = ch_model    // path: scVI model directory
    ref_paths  = ch_refs     // channel: [ ref.h5ad, ... ]
    versions   = ch_versions // channel: [ versions.yml ]
}
