#!/bin/bash

nextflow sc-annotate.nf \
	-profile conda \
	-params-file params.hs.json \
	--study_names tests/study_names_human.txt \
	-process.executor local \
	--process_samples false \
	--mask false \
	-resume
	
