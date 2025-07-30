include { DOWNLOAD_STUDIES } from "${projectDir}/modules/processes/download_studies.nf"

workflow DOWNLOAD_STUDIES_SUBWF {
    
    take:
    // One input must be null
    // Determines whether to download studies or use pre-existing path
    study_names
    study_paths

    main:
    if (study_names) {
        // change to accept command line input or params input separated by space
        
	// Define the study name
    Channel
        .from(params.study_names.split(/\s+/))
        .set { study_names }
            
            //.fromPath(params.study_names)
            //.flatMap { file -> file.readLines().collect { it.trim() } }
            //.set { study_names }

        DOWNLOAD_STUDIES(study_names)
            .set { study_channel }

    } else if (study_paths) {
        // change to accept command line input or params input
        // separated by space

        def paths = params.study_paths

        // Support both repeated flags and space‑separated string
        if (!(paths instanceof List)) {
            paths = paths.toString().split(/\s+/)
        }


        Channel
            .from( paths.collect { file(it) } ) // convert to Path objects
            .set { study_paths }

        study_paths.view()
        // get study names from each path
        study_paths.map { path ->
            def name = path.getName()
            [name, path]
        }.set { study_channel }
            //.fromPath(params.study_paths)
            //.flatMap { path ->
                //def results = []
                //path.eachDir { dir -> results << [dir.name, dir] }
                //return results
            //}
    } else {
        exit 1, "Error: You must provide either 'study_names' or 'study_paths'."
    }

    emit: study_channel
    
}
