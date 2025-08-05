include { DOWNLOAD_STUDIES } from "${projectDir}/modules/processes/download_studies.nf"

workflow DOWNLOAD_STUDIES_SUBWF {

    take:
    study_names
    study_paths

    main:

    if (study_names) {
        // read from file or space‑separated list
        study_channel = (
            file(study_names).exists() && file(study_names).isFile()
                ? Channel.from(
                      file(study_names)
                          .readLines()
                          .collect { it.trim() }
                          .findAll { it }
                          .collect { it }
                  )
                : Channel.from(
                      study_names
                          .toString()
                          .split(/\s+/)
                          .findAll { it }
                          .collect { it }
                  )
        )

        study_channel = DOWNLOAD_STUDIES(study_channel)

    }
    else if (study_paths) {
        // read from file or space‑separated list
        study_channel = (
            file(study_paths).exists() && file(study_paths).isFile()
                ? Channel.from(
                      file(study_paths)
                          .readLines()
                          .collect { it.trim() }
                          .findAll { it }
                          .collect { path ->
                              def p = file(path)
                              [p.getName(), p]
                          }
                  )
                : Channel.from(
                      study_paths
                          .toString()
                          .split(/\s+/)
                          .findAll { it }
                          .collect { path ->
                              def p = file(path)
                              [p.getName(), p]
                          }
                  )
        )
    }
    else {
        exit 1, "Error: You must provide either 'study_names' or 'study_paths'."
    }

    emit:
    study_channel
}
