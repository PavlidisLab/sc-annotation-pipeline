#!/bin/bash

nextflow sc-annotate.nf \
	-profile conda \
	-params-file params.mm.json \
	--study_names GSE152715.1 \
	-process.executor local \
	--process_samples true \
	-resume \
	--mask false 
	
