#!/bin/bash

nextflow main.nf \
	-profile conda \
	-params-file params.hs.json \
	--study_names /space/grp/rschwartz/rschwartz/cell_annotation_cortex.nf/tests/study_names_human.txt \
	-process.executor local \
	--process_samples false \
	--mask false \
	-resume
	
