#!/bin/bash
# Amy Campbell
# Staphyloxanthin paper analyses 
#
# Step 1:Trimming nextera transposase adapters & low quality bases off
#  the raw Illumina reads
# outputs trimmed reads into WholeGenomeSequencing/TrimmedReads/

conda activate TrimmingEnvironment

mkdir -p TrimmedReads

for readsfile1 in ../Data/WGSData/RawIllumina/*R1 ; do
	basestring=$(basename $readsfile1)
	r1ext="_R1.fastq"
	r2ext="_R2.fastq"
	blankstring=""

	prefixUsed=${basestring/$r1ext/$blankstring} 
	readsfile2=${readsfile1/$r1ext/$r2ext}
	
	trimmedext="trimmed"
	basefilename=$prefixUsed

	trim_galore --paired --nextera --clip_R1 10 --clip_R2 10 --length 70 --stringency 10 --basename $prefixUsed --output_dir TrimmedReads $readsfile1 $readsfile2
done
