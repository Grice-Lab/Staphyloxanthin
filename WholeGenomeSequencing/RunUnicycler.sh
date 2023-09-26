#!/bin/bash
# Amy Campbell
# Staphyloxanthin paper analyses 
# Unicycler calls for both short read only and hybrid assembly
# requires trimming steps to be performed before use

relativepathData="../Data/WGSData/"
fullpathData=$(echo "$(cd "$(dirname "$relativepathData")"; pwd)/$(basename "$relativepathData")")


InputTrimmedNanopore=$fullpathData"/TrimmedONP/"
InputTrimmedIllumina=$fullpathData"/TrimmedReads/"

OutputIlluminaOnly=$fullpathData"/IlluminaAssemblies/"
OutputHybrid=$fullpathData"/HybridAssemblies/"

OutputIlluminaContigs=$fullpathData"/IlluminaAssemblies/Contigs/"
OutputHybridContigs=$fullpathData"HybridAssemblies/Contigs/"

mkdir -p $OutputIlluminaContigs
mkdir -p $OutputIlluminaOnly
mkdir -p $OutputHybrid
mkdir -p $OutputHybridContigs

conda activate EAGenomeEnv

# Illumina-only assemblies
##########################
# actually change to reflect trimmed names
for FwdReads in $InputTrimmedIllumina*trimmed_val_1.fastq ; do
	revsuffix="trimmed_val_2.fastq"
	fwdsuffix="trimmed_val_1.fastq"
	blank=""

	RevReads=${FwdReads/$fwdsuffix/$revsuffix}
	
	baseString=$(basename $FwdReads)
	JustID=${baseString/$fwdsuffix/$blank}
	
	unicycler -1 $FwdReads -2 $RevReads -o $OutputIlluminaOnly$JustID
	cp $OutputIlluminaOnly$JustID"/assembly.fasta" $OutputIlluminaContigs$JustID".fasta"
done

for reads in $InputTrimmedONP*_trimmedONP.fastq ; do
	baseString=$(basename $reads)
	replacestring="_trimmedONP.fastq"
	blank=""
   	JustID=${baseString/$replacestring/$blank}
	fwdextShort="trimmed_val_1.fastq"
	revextShort="trimmed_val_2.fastq"

	fwdshort=$InputTrimmedIllumina$JustID$fwdextShort
	revshort=$InputTrimmedIllumina$JustID$revextShort


	unicycler -1 $fwdshort -2 $revshort -l $reads -o $OutputHybrid$IDonly
	cp $OutputHybrid$IDonly"/assembly.fasta" $OutputHybridContigs$JustID".fasta"
done


