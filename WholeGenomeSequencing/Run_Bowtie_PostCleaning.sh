#!/bin/bash
# Amy Campbell
# Adapted 9/22 for Staphyloxanthin paper
# Shell script to map each genome's trimmed raw sequencing reads 
# Back to the assembled + cleaned contigs for that genome
# for use downstream in TestAdmixtures.sh

relativepathData="../Data/"
fullpathData=$(echo "$(cd "$(dirname "$relativepathData")"; pwd)/$(basename "$relativepathData")")
fullpathBowtie=$fullpathData"/BowtieDatabases/"
conda activate BowtieEnv


export BOWTIE2_INDEXES=$fullpathBowtie


outputfolder=$fullpathData"/WGSData/Coverages_Cleaned/"

mkdir -p $outputfolder

Raw_FastQs=$fullpathData"/WGSData/TrimmedReads/"
FinalContigs=$fullpathData"/WGSData/FinalIlluminaContigs/"

for filename in $FinalContigs*.fasta; do 
	filenamestring=$(basename $filename)
	ext="_cleaned.fasta"
	blank=""
	noext=${filenamestring/$ext/$blank}
	
	bowtie2-build $filename $noext
	mv *.bt2 $fullpathBowtie

	sam_ext=".sam"
	fwd_Ext="trimmed_val_1.fastq"
	rev_Ext="trimmed_val_2.fastq"


	fwd_reads=$Raw_FastQs$noext$fwd_Ext
	rev_reads=$Raw_FastQs$noext$rev_Ext

	samfile=$noext$sam_ext
	
	bowtie2 -x $noext -1 $fwd_reads -2 $rev_reads -S $outputfolder$samfile

done
