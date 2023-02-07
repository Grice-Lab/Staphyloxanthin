# Merge demultiplexed fwd and reverse reads across lanes within replicates

demuxedreadsprefix=/home/acampbe/DFU/data/RNASeq/RawReadsDemux/
outputpath="/home/acampbe/DFU/data/RNASeq/MergedLanes/"

mkdir -p $outputpath

for prefix in ARM1 ARM2 ARM3 ARM4 ARM5 ARM7 ARM8 ARM9 ARM10 ARM11 ARM12; do
 fwdreadsext="_R1.fastq"
 reversereadsext="_R2.fastq"
 
 outputfwd=$outputpath$prefix$fwdreadsext
 outputrev=$outputpath$prefix$reversereadsext

 touch $outputfwd
 touch $outputrev

 for samplefile in $demuxedreadsprefix$prefix* ; do
 
    forwardreads=$samplefile/*R1_001.fastq.gz
    revreads=$samplefile/*R2_001.fastq.gz
	
    gunzip $forwardreads
    gunzip $revreads
  

    cat $samplefile/*R1_001.fastq >> $outputfwd
    cat $samplefile/*R2_001.fastq >> $outputrev 

    
  
 done 


done
