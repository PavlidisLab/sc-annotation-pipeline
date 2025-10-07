#!/bin/bash

nextflow sc-annotate.nf \
	-profile conda \
	-params-file params.mm.json \
	--study_names "GSE124952" \
	-process.executor local \
	--process_samples \
	--use_staging \
	-resume \
	--mask 
	
