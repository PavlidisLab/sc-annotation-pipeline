#!/bin/bash

nextflow sc-annotate.nf \
	-profile conda \
	-params-file params.mm.json \
	--study_paths "GSE223423 GSE295078" \
	-process.executor local \
	--process_samples false \
	-resume \
	--mask false 
	
