#!/bin/bash

nextflow sc-annotate.nf \
	-profile conda \
	-params-file params.mm.json \
	--study_paths /space/grp/rschwartz/rschwartz/cell_annotation_cortex.nf/tests/study_paths_mouse.txt \
	-process.executor local \
	--process_samples false \
	-resume \
	--mask false 
	
