
process PROCESS_QUERY_SAMPLE {

	tag "$study_name"

	// This process is used to process query datasets using a pre-trained model.
	// It takes a study name and path, processes the data, and outputs processed and raw data files.
	// The model path is provided as an input parameter.


    input:
    val model_path
    tuple val(study_name), val(study_path)

    output:
    tuple val("${study_name}"), path("${study_name}.h5ad"), emit: processed_query
    tuple val("${study_name}"), path("${study_name}_raw.h5ad"), emit: raw_query

        
    script:


    """

    python $projectDir/bin/process_query_samples.py \\
                            --model_path ${model_path} \\
                            --query_name ${study_name} \\
                            --query_path ${study_path} \\
                            --seed ${params.seed}
    """

}