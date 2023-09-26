#!/bin/bash
# Amy Campbell
# Modified for Staphyloxanthin paper 
# Using BCFtools and samtools to call SNPs
# Based on Raven et al. 2020 paper 

conda activate BowtieEnv

relativepathData="../Data/"
fullpathData=$(echo "$(cd "$(dirname "$relativepathData")"; pwd)/$(basename "$relativepathData")")

OutputPath=$fullpathData"/WGSData/BCF_Variant_Calling/" 
PathToSams=$fullpathData"/WGSData/Coverages_Cleaned/"


mkdir -p $OutputPath

#samtools mpileup -uf tmp/reference.fa input/alignments/*.bam > tmp/variants/raw_calls.bcf

for filename in /project/grice/storage/DFUShortReads2022/FinalContigs/*.fasta; do 
        filenamestring=$(basename $filename)
        ext="_cleaned.fasta"
        blank=""
        noext=${filenamestring/$ext/$blank}
	bamext=".bam"
	samext=".sam"
	baiext=".bai"
	rawvcfext=".raw.vcf"
	bcfext=".bcf"
	MAFext="_MAF.bcf"

	# Samtools index	
	################
	samtools faidx $filename
	
	# Convert to bam 
	#################
	samtools view -Sb $PathToSams$noext$samext | samtools sort - > $PathToSams$noext$bamext
	samtools index $PathToSams$noext$bamext
	 
	#  
	################ 
	samtools mpileup --skip-indels -uf $filename $PathToSams$noext$bamext > $PathToSams$noext$rawvcfext
	bcftools call $PathToSams$noext$rawvcfext --ploidy 1 --skip-variants indels -mv > $OutputPath$noext$bcfext
	bcftools view -i 'MAF > 0.05' $OutputPath$noext$bcfext -Ov -o  $OutputPath$noext$MAFext

done	
	
