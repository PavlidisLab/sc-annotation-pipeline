#!/bin/bash

nextflow sc-annotate.nf \
	-profile conda \
	-params-file params.mm.json \
	--study_paths /space/grp/rschwartz/rschwartz/get_gemma_data.nf/study_names_mouse.txt_author_true_process_samples_true/mex/GSE247339.2 \
	-process.executor local \
	--process_samples false \
	-resume \
	--mask false 
	
