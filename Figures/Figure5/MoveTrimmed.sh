#/home/acampbe/DFU/data/RNASeq
#ARM9_val_2.fq


for trimmedfile in /home/acampbe/DFU/data/RNASeq/Trimmed/*_val_1*; do
	
	fwdext_current="_val_1.fq"
        revext_current="_val_2.fq"
        fwdext_new="trimmed_R1.fastq"
        revext_new="trimmed_R2.fastq"

	currentfilefwd=$trimmedfile
	currentfilerev=${trimmedfile/$fwdext_current/$revext_current}

	
	newfwd=${trimmedfile/$fwdext_current/$fwdext_new}
	newrev=${trimmedfile/$fwdext_current/$revext_new}

	mv $currentfilefwd $newfwd
	mv $currentfilerev $newrev
	
done
