include { DOWNLOAD_STUDIES } from "${projectDir}/modules/processes/download_studies.nf"

workflow DOWNLOAD_STUDIES_SUBWF {
    
    take:
    // One input must be null
    // Determines whether to download studies or use pre-existing path
    study_names
    studies_path

    main:
    if (study_names) {
        Channel
            .fromPath(params.study_names)
            .flatMap { file -> file.readLines().collect { it.trim() } }
            .set { study_names }

        DOWNLOAD_STUDIES(study_names)
            .set { study_channel }

    } else if (studies_path) {
        study_channel = Channel
            .fromPath(params.studies_path)
            .flatMap { path ->
                def results = []
                path.eachDir { dir -> results << [dir.name, dir] }
                return results
            }
    } else {
        exit 1, "Error: You must provide either 'study_names' or 'studies_path'."
    }

    emit: study_channel
    
}
