#!/bin/bash

nextflow main.nf \
	-profile conda \
	-params-file params.mm.json \
	--study_names tests/study_names_mouse.txt \
	-process.executor local \
	-resume 
	
