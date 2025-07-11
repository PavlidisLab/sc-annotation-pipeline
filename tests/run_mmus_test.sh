#!/bin/bash

nextflow sc-annotate.nf \
	-profile conda \
	-params-file params.mm.json \
	--study_names /space/grp/rschwartz/rschwartz/cell_annotation_cortex.nf/study_names_mouse.txt \
	-process.executor slurm \
	--process_samples true \
	-resume
	
