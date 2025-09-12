#!/bin/bash

nextflow main.nf \
	-profile conda \
	-params-file params.hs.json \
	--study_names tests/study_names_human.txt \
	-process.executor local \
	-resume
	
