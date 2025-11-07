#!/bin/bash

nextflow sc-annotate.nf \
	-profile conda \
	-params-file params.hs.json \
	--study_paths "/space/grp/rschwartz/rschwartz/cell_annotation_cortex.nf/results/homo_sapiens_subsample_ref_500_2025-09-26_13-37-29/mex/GSE278619" \
	-process.executor local \
	--process_samples false \
	--mask false \
	-resume
	
