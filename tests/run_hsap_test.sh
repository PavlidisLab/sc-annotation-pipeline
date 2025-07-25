#!/bin/bash

nextflow sc-annotate.nf \
	-profile conda \
	-params-file params.hs.json \
	--studies_path /space/grp/rschwartz/rschwartz/cell_annotation_cortex.nf/results/homo_sapiens_subsample_ref_500_2025-07-21_13-55-05/mex/Velmeshev_et_al.2 \
	-process.executor local \
	--process_samples false \
	--mask false \
	-resume
	
