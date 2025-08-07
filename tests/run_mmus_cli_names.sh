#!/bin/bash

nextflow sc-annotate.nf \
	-profile conda \
	-params-file params.mm.json \
	--study_names "GSE124952 GSE185454" \
	-process.executor local \
	--process_samples false \
	-resume \
	--mask 
	
