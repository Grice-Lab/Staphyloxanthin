#!/bin/bash
# Jan 22 
# Copied from general CoverageAnalyses script called Run_Bowtie_DORNs.sh 
# Shell script to map each DORN genome's raw sequencing reads 
# Back to the assembled contigs for that genome 
# In order to assess breadth of coverage, depth
# This is used for the first round of genome cleaning

source /home/acampbe/software/miniconda3/bin/activate BowtieEnv

mkdir -p /home/acampbe/DFU/data/AlignmentCoverage/CompleteCoverageData/BowtieDatabases/

export BOWTIE2_INDEXES=/home/acampbe/DFU/data/AlignmentCoverage/CompleteCoverageData/BowtieDatabases/


# Build bowtie2 databases for every genome in the assembly folder
outputfolder="/home/acampbe/DFU/data/AlignmentCoverage/CompleteCoverageData/InitialAlignment/"

mkdir -p $outputfolder


Trimmed_FastQs="/project/grice/storage/DFUShortReads2022/trimmedreads/"
Assemblies="/project/grice/storage/DFUShortReads2022/assemblies"

for filename in $Assemblies/*.fasta; do 
	filenamestring=$(basename $filename)
	ext="_contigs.fasta"
	blank=""
	noext=${filenamestring/$ext/$blank}
	
	bowtie2-build $filename $noext
	mv *.bt2 /home/acampbe/DFU/data/AlignmentCoverage/CompleteCoverageData/BowtieDatabases/
	
	sam_ext=".sam"
	fwd_Ext="trimmedgalore_val_1.fastq"
	rev_Ext="trimmedgalore_val_2.fastq"


	fwd_reads=$Trimmed_FastQs$noext$fwd_Ext
	rev_reads=$Trimmed_FastQs$noext$rev_Ext

	samfile=$noext$sam_ext
	
	bowtie2 -x $noext -1 $fwd_reads -2 $rev_reads -S $outputfolder$samfile

done
