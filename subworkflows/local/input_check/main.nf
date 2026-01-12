/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS AND PREPARE STUDY CHANNELS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { DOWNLOAD_STUDIES } from '../../../modules/local/download_studies/main'

workflow INPUT_CHECK {

    take:
    input         // path: samplesheet CSV file (optional)
    study_names   // string: space-separated list or file path (optional)
    study_paths   // string: space-separated list or file path (optional)
    use_staging   // boolean: use Gemma staging environment

    main:
    ch_versions = Channel.empty()

    //
    // Option 1: Samplesheet input (--input)
    //
    if (input) {
        Channel
            .fromPath(input, checkIfExists: true)
            .splitCsv(header: true, strip: true)
            .map { row ->
                // Validate row has sample column
                if (!row.sample) {
                    error "Error: Samplesheet row missing required 'sample' column: ${row}"
                }
                // Return tuple with sample, study_name, study_path
                def sample = row.sample.trim()
                def study_name = row.study_name ? row.study_name.trim() : null
                def study_path = row.study_path ? row.study_path.trim() : null

                if (!study_name && !study_path) {
                    error "Error: Row '${sample}' must have either 'study_name' or 'study_path'"
                }

                [sample, study_name, study_path]
            }
            .branch { sample, study_name, study_path ->
                download: study_name && !study_path
                    return study_name
                local: study_path
                    return [sample, file(study_path)]
            }
            .set { ch_samplesheet }

        // Download studies that need downloading
        DOWNLOAD_STUDIES(
            ch_samplesheet.download,
            use_staging
        )

        // Combine downloaded and local studies
        ch_studies = ch_samplesheet.local.mix(DOWNLOAD_STUDIES.out.study_dir)
    }
    //
    // Option 2: Study names (--study_names) - download from Gemma
    //
    else if (study_names) {
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
    //
    // Option 3: Study paths (--study_paths) - local directories
    //
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
    //
    // No input provided
    //
    else {
        error "Error: You must provide one of: '--input' (samplesheet), '--study_names', or '--study_paths'."
    }

    emit:
    studies  = ch_studies   // channel: [ sample/study_name, study_path ]
    versions = ch_versions  // channel: [ versions.yml ]
}
