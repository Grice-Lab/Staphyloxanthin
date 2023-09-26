#!/bin/bash
# Amy Campbell
# Staphyloxanthin paper 2022
# Combining the output of coverage_stats.py from each of the Illumina genome assemblies (2019 & 2022)


relativepathData="../Data/"
fullpathData=$(echo "$(cd "$(dirname "$relativepathData")"; pwd)/$(basename "$relativepathData")")

inputTSVs=$fullpathData"/WGSData/CoverageStatsContigTSVs/"
outputfile=$fullpathData"/WGSData/AllCoverageStats.tsv"

echo "identifer    length    mapped    placed    min_cov    max_cov    mean_cov" > $outputfile


for tsvfile in $inputTSVs*.tsv ; do

        justfname=$(basename $tsvfile)
        ext=".tsv"
        blank=""
        underscore="_"
        basenamefile=${justfname/$ext/$blank}
        addition=$basenamefile$underscore

        # Add genome identifier to the beginning of each line to delineate contigs by genome of origin
	
	awk -v prefix="$addition" '{print prefix $0}' $tsvfile > temp.csv

        # delete header from this temporary file so that we can assign a single header (above) to the concatenated file
        sed -i '1d' temp.csv 

        # add to output file
	cat temp.csv >> $outputfile
	rm temp.csv
done

