#!/bin/bash

nextflow sc-annotate.nf \
	-profile conda \
	-params-file params.hs.json \
	--study_paths /space/grp/rschwartz/rschwartz/get_gemma_data.nf/psychEncode.txt_author_true_sample_split_true/mex_dirs/Velmeshev_et_al.1 \
	-process.executor local \
	--process_samples false \
	--mask false \
	-resume
	
