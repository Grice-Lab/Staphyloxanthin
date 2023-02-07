# Amy Campbell 
# 2022
# Step 1 RNASeq on Two SA strains: FastQC

source /home/acampbe/software/miniconda3/bin/activate TrimmingEnvironment

for filename in $READSPATH/*_R1.fastq.gz; do
        run1=$filename  
        run2=${run1/$r1string/$r2string}
        
        fastqc -o $OUTPUT $run1 -f fastq
        fastqc -o $OUTPUT $run2 -f fastq

done    

