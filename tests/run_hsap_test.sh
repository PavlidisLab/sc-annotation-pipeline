#!/bin/bash

nextflow sc-annotate.nf \
	-profile conda \
	-params-file params.hs.json \
	--study_names /space/grp/rschwartz/rschwartz/cell_annotation_cortex.nf/study_names_human.txt \
	-process.executor slurm \
	-resume
	
