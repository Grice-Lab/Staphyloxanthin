#!/bin/bash
# Amy Campbell
# Staphyloxanthin paper analyses
#
# Step 3: run porechopper to trim adapters from the ONP reads, quality-filter

CONDA_ACTIVATE_PATH=$1
READSPATH=$2
OUTPUT=$3

conda activate  ONPEnv

mkdir -p $OUTPUT

outputextension="_trimmedONP.fastq"

for ONPreads in 
	
	porechop -i $READSPATH$FNAME -o $OUTPUT$ID$outputextension

done < $MAPPING
