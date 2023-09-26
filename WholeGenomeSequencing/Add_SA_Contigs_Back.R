# Amy Campbell
# Jan 2022
# As part of the 'contig cleaning' portion of my sequence data processing 
# pipeline for this project, this script reads the tabular output of the blast search 
# for each 'questionable' contig (either a high-coverage or low-coverage outlier for its
# genome of origin)

# Each time this script is called, it reads in the file located at the path 
# input with the call
# 

# INPUTS
#########
# Tabular blast output in the form GenomeName_FollowUpContigs.tab 
# I use the tabular output from the blast searches I performed on 
# the 'follow up contigs' identified as outliers by depth

# EXAMPLE CALL
#############
# Rscript Add_SA_contigs_Back.R /project/grice/storage/HiSeq/WGS/HiSeq_19/AssemblyFastas/DFU100_followupContigs/BlastOutput/DORN1447_FollowUpContigs.tab

library("dplyr")
library("stringr")

args <- commandArgs(trailingOnly=T)
filepath=args[1]

#filepath="/Users/amycampbell/DORN1447_FollowUpContigs.tab"
filename = basename(filepath)
genomename = stringr::str_remove(filename, "_FollowUpContigs.tab")

blastoutput = read.csv(filepath, sep='\t',col.names = c("contig_num","taxIDs", "Species", "accessions", "pctID", "length", "gapopen", "qstart", "qend", "eval","bitscore"))

# Highest bit score hit for each followed-up contig in this genome. 
TopScoringHits = blastoutput %>% group_by(contig_num) %>% slice(which.max(bitscore))

# Sets a logical var to true if the species of the top hit for a contig has 'Staphylococcus'
# Choosing to keep 'staphylococcus' and not aureus in particular because I want to allow for highly conserved 
# staphylococcal regions, staph phages 
TopScoringHits = TopScoringHits %>% rowwise() %>% mutate(Keep = (grepl("Staphylococcus", Species) || grepl("staphylococcal", Species)))

# As long as none are 
NonStaphs = TopScoringHits %>% filter(!Keep)

# If there are any contigs that needn't be added back 
if(dim(NonStaphs)[1]!=0){
  write.csv(NonStaphs[c("contig_num", "Species")], file=paste0(genomename, "_RemovedAfterBlast.txt"))
}else{
  print(paste0("None_", genomename))
}


