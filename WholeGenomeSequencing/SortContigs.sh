#!/bin/bash
# amy campbell 
# shellscript to run SortingContigsFate.py
# Modified for September 2022 Staphyloxanthin paper

conda activate BowtieEnv

relativepathData="../Data/"
fullpathData=$(echo "$(cd "$(dirname "$relativepathData")"; pwd)/$(basename "$relativepathData")")
 

assemblies=$fullpathData"/WGSData/IlluminaAssemblies/Contigs/"
fates=$fullpathData"/WGSData/ContigFates.tsv"
clean=$fullpathData"/WGSData/ShortAssemblyCleanContigs/"
follow=$fullpathData"/WGSData/ShortAssemblyFollowUpContigs/"

mkdir -p $follow
mkdir -p $clean


python SortingContigsFate.py $assemblies $fates $clean $follow

