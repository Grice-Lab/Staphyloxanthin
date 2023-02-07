# Amy Campbell
# Align RNAseq paired end, trimmed reads (ARM1..12) to reference
# 04/2022

source /home/acampbe/software/miniconda3/bin/activate BowtieEnv


bowtiepath="/home/acampbe/DFU/data/RNASeq/BTIndices/"
readspath="/home/acampbe/DFU/data/RNASeq/Trimmed"
outputpath="/home/acampbe/DFU/data/RNASeq/BTAlign/"
refpath="/home/acampbe/DFU/data/RNASeq/DORN925Reference.fasta"

mkdir -p $bowtiepath
export BOWTIE2_INDEXES=$bowtiepath

# Build bowtie index
####################
refname="DORN925_Ref"

bowtie2-build $refpath $refname
mv *.bt2 $bowtiepath

forwardext="trimmed_R1.fastq"
revext="trimmed_R2.fastq"
noext=""
samext=".sam"
bamext=".bam"
sortedbamext="_sorted.bam"

for filename in $readspath/*R1.fastq ; do 
	base=$(basename $filename)
	prefixname=${base/$forwardext/$noext}
	
	forwardfile=$filename
	reversefile=${filename/$forwardext/$revext}

	bowtie2 -x $refname -1 $forwardfile -2 $reversefile -S $outputpath$prefixname$samext
	samtools view -bS $outputpath$prefixname$samext > $outputpath$prefixname$bamext	
	samtools sort $outputpath$prefixname$bamext > $outputpath$prefixname$sortedbamext
done

