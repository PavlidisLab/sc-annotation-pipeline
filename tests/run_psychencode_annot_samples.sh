#!/bin/bash

nextflow sc-annotate.nf \
	-profile conda \
	-params-file params.hs.json \
	--study_paths psychencode_paths_large.txt \
	-process.executor slurm \
	--process_samples true \
	--mask false \
	-work-dir work \
	--use_staging true \
	-resume
	
