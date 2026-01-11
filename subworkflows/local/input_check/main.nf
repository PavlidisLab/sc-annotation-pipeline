/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS AND PREPARE STUDY CHANNELS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { DOWNLOAD_STUDIES } from '../../../modules/local/download_studies/main'

workflow INPUT_CHECK {

    take:
    study_names   // string: space-separated list or file path
    study_paths   // string: space-separated list or file path
    use_staging   // boolean: use Gemma staging environment

    main:
    ch_versions = Channel.empty()

    if (study_names) {
        // Read from file or space-separated list
        study_channel = (
            file(study_names).exists() && file(study_names).isFile()
                ? Channel.from(
                      file(study_names)
                          .readLines()
                          .collect { it.trim() }
                          .findAll { it }
                  )
                : Channel.from(
                      study_names
                          .toString()
                          .split(/\s+/)
                          .findAll { it }
                  )
        )

        // Download studies from Gemma
        DOWNLOAD_STUDIES(study_channel, use_staging)
        ch_studies = DOWNLOAD_STUDIES.out.study_dir
    }
    else if (study_paths) {
        // Read from file or space-separated list
        ch_studies = (
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
        error "Error: You must provide either 'study_names' or 'study_paths'."
    }

    emit:
    studies  = ch_studies   // channel: [ study_name, study_path ]
    versions = ch_versions  // channel: [ versions.yml ]
}
