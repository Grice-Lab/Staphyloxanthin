#!/bin/bash
# Amy Campbell
# Staphyloxanthin paper 2022
# Do a blast search of each contig that was selected for 'followup' based on its being a depth of coverage outlier or a length outlier

# activate condas environment with blast+ installed
conda activate BlastEnv

relativepathData="../Data/"
fullpathData=$(echo "$(cd "$(dirname "$relativepathData")"; pwd)/$(basename "$relativepathData")")


# Set path to  NT database from ncbi (version I used  was updated 6/22/2020)
# SET THIS(/home/acampbe/DownloadedDatabases/BlastDBs/ for me):
blastpath=""

export BLASTDB=$blastpath

followup=$fullpathData"/WGSData/ShortAssemblyFollowUpContigs/"

outputfolder=$fulldata"/WGSData/BlastOutput/"
mkdir -p $outputfolder


for f in $followup*.fasta ; do
	ext=".fasta"
	blank=""
	outext=".tab"
	justfname=$(basename $f)
	noext=${justfname/$ext/$blank}

		
	blastn -query $f -db "nt" -outfmt "6 qseqid staxids sscinames sseqid pident length gapopen qstart qend evalue bitscore" -out $outputfolder$noext$outext

done

