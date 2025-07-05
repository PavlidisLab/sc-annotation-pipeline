include { PROCESS_QUERY_SAMPLE } from "${projectDir}/modules/processes/process_query_samples.nf"
include { PROCESS_QUERY_COMBINED } from "${projectDir}/modules/processes/process_query_combined.nf"

workflow PROCESS_QUERY_SUBWF {
    take:
	study_channel
	model_path
	

	main:
	// Process each query by relabeling, subsampling, and passing through scvi model

	if params.process_by_sample {
		// split study_channel into individual samples
		study_channel.flatMap { study_name, study_dir ->
		    def results = []
		    study_dir.eachDir { dir -> 
				results << [dir.name, dir] 
			}
		    return results
		}.set { study_channel }

	    // If process_by_sample is true, we will
	    // process each query sample separately and use a different process
	    PROCESS_QUERY_SAMPLE(model_path, study_channel)
	} else {
	    // Process each query without subsampling
	    PROCESS_QUERY_COMBINED(model_path, study_channel)
	}

}
