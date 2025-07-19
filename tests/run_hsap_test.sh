#!/bin/bash

nextflow sc-annotate.nf \
	-profile conda \
	-params-file params.hs.json \
	--studies_path //home/bxu/sc-annotation-pipeline/work/74/0e54911c90f1ab60613cd463f5d73c \
	-process.executor slurm \
	--process_samples true \
	--mask false \
	-resume
	
