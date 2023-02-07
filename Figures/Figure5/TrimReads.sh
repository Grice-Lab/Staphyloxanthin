# Amy Campbell 
# 2022
# Step 1 RNASeq on Two SA strains: FastQC

source /home/acampbe/software/miniconda3/bin/activate TrimmingEnvironment

ReadsPath="/home/acampbe/DFU/data/RNASeq/MergedLanes/"
Output="/home/acampbe/DFU/data/RNASeq/Trimmed"
mkdir -p $Output

for filename in $ReadsPath/*_R1.fastq; do
	r1string="_R1.fastq"
	r2string="_R2.fastq"
	blankstring=""

        run1=$filename 
        run2=${run1/$r1string/$r2string}
        basefname=$(basename $filename)
	basefnameNoExt=${basefname/$r1string/$blankstring}

        trim_galore --paired --clip_R1 1 --clip_R2 1 --basename $basefnameNoExt --output_dir $Output $run1 $run2

done    

