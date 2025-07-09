include { PROCESS_QUERY_SAMPLE } from "${projectDir}/modules/processes/process_query_samples.nf"
include { PROCESS_QUERY_COMBINED } from "${projectDir}/modules/processes/process_query_combined.nf"

workflow PROCESS_QUERY_SUBWF {

    take:
	study_channel
	model_path
	process_samples

	main:
		def processed_queries
		def raw_queries
		if (process_samples == true) {
			// Split study_channel into individual samples
			expanded_channel = study_channel.flatMap { study_name, study_dir ->
					def results = []
					study_dir.eachDir { dir -> results << [dir.name, dir] }
					return results
				}
			// Process each query sample separately
			PROCESS_QUERY_SAMPLE(model_path, expanded_channel)
			raw_queries = PROCESS_QUERY_SAMPLE.out.raw_query
			processed_queries = PROCESS_QUERY_SAMPLE.out.processed_query
		} else  if (process_samples == false) {
			// Process each query without subsampling
			PROCESS_QUERY_COMBINED(model_path, study_channel)
			raw_queries = PROCESS_QUERY_COMBINED.out.raw_query
			processed_queries = PROCESS_QUERY_COMBINED.out.processed_query
			
		} else {
			error "Invalid value for process_samples: ${process_samples}. Expected true or false."
		}
	emit:
		processed_queries  // Processed queries after subsampling and relabeling
		raw_queries
}