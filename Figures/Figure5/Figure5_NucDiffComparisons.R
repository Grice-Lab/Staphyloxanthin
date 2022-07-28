# Amy Campbell
# 2022
# Programmatic comparison of NucDiff results for DORN925 vs. each of the other 13 isolates in Patient # 141
# Previously done by looking for presence/absence of each variant found in DORN1000 compared to DORN925
# (where DORN925 was 'reference' and DORN1000 was 'query')
# Now, doing with DORN1088 because it's the one we used in the mouse experiments, RNAseq

library(dplyr)

# Get list of variants & locations in DORN1088 compared to DORN925
##################################################################

SNPs1088 = read.csv2("data/nucdiff925/NucDiff_925_1088/results/DORN925_DORN1088Comparison_ref_snps.gff",comment.char = "#", sep="\t",header=F)

colnames(SNPs1088) = c("Contig", "NucDiffVersion","SO", "Start", "End", "Blank1", "Blank2", "Blank3", "Description")
SNPs1088 = SNPs1088 %>% select(Start, End, Description)


Structural1088 = read.csv2("data/nucdiff925/NucDiff_925_1088/results/DORN925_DORN1088Comparison_ref_struct.gff",comment.char = "#", sep="\t",header=F)
colnames(Structural1088) = c("Contig", "NucDiffVersion","SO", "Start", "End", "Blank1", "Blank2", "Blank3", "Description")
Structural1088 = Structural1088 %>% select(Start, End, Description)

Additional1088 = read.csv2("data/nucdiff925/NucDiff_925_1088/results/DORN925_DORN1088Comparison_ref_additional.gff",comment.char = "#", sep="\t",header=F)
colnames(Additional1088) = c("Contig", "NucDiffVersion","SO", "Start", "End", "Blank1", "Blank2", "Blank3", "Description")
Additional1088 = Additional1088 %>% select(Start, End, Description)

Variants1088 = rbind(SNPs1088,Structural1088)
Variants1088= rbind(Variants1088, Additional1088)

genomevector=c("DORN1088", "DORN1000", "DORN1038", "DORN1194", "DORN1061", "DORN1037", "DORN999", "DORN976", "DORN933", "DORN929", "DORN882", "DORN881", "DORN880")
Results=Variants1088
for(GenomeString in genomevector){
  ComparisonString=paste0("DORN925_", GenomeString)
  JustNumbers = stringr::str_remove_all(ComparisonString, "DORN")
  print(ComparisonString)
  SNPFile =paste("data/nucdiff925/NucDiff_", JustNumbers, "/results/", ComparisonString, "Comparison_ref_snps.gff",sep="")
  AdditionalFile =paste("data/nucdiff925/NucDiff_", JustNumbers, "/results/", ComparisonString, "Comparison_ref_additional.gff",sep="")
  StructFile =paste("data/nucdiff925/NucDiff_", JustNumbers, "/results/", ComparisonString, "Comparison_ref_struct.gff",sep="")
  
  #SNPdf = read.csv2(SNPFile,comment.char = "#", sep="\t",header=F)
  SNPdf = tryCatch({read.csv2(SNPFile,comment.char = "#", sep="\t",header=F)}, error=function(errmessage){
    message("No lines in that file")
    return(data.frame())
  })
  if(nrow(SNPdf) > 0){
  colnames(SNPdf) = c("Contig", "NucDiffVersion","SO", "Start", "End", "Blank1", "Blank2", "Blank3", "Description")
  SNPdf = SNPdf %>% select(Start, End, Description)}
  
  #structDF = read.csv2(StructFile,comment.char = "#", sep="\t",header=F)
  structDF = tryCatch({read.csv2(StructFile,comment.char = "#", sep="\t",header=F)}, error=function(errmessage){
    message("No lines in that file")
    return(data.frame())
  })
  
  if(nrow(structDF)>0){
    colnames(structDF) = c("Contig", "NucDiffVersion","SO", "Start", "End", "Blank1", "Blank2", "Blank3", "Description")
    structDF = structDF %>% select(Start, End, Description)
  }

  #AdditionalDF = read.csv2(AdditionalFile,comment.char = "#", sep="\t",header=F)
  AdditionalDF = tryCatch({read.csv2(AdditionalFile,comment.char = "#", sep="\t",header=F)}, error=function(errmessage){
    message("No lines in that file")
    return(data.frame())
  })
  if(nrow(AdditionalDF)>0){
    colnames(AdditionalDF) = c("Contig", "NucDiffVersion","SO", "Start", "End", "Blank1", "Blank2", "Blank3", "Description")
    AdditionalDF = AdditionalDF %>% select(Start, End, Description)
  }

  
  VariantDF = rbind(SNPdf, AdditionalDF)
  VariantDF = rbind(VariantDF, structDF)
  
  VariantVector = rep(NA, nrow(Variants1088))
  for(rownum in 1:nrow(Variants1088)){
    row=Variants1088[rownum,]
    startpresent = (row["Start"] %in% VariantDF$Start)
    endpresent = (row["End"] %in% VariantDF$End)
    if((startpresent & endpresent )){
      VariantVector[rownum] <- 1
    }else if(startpresent | endpresent){
      VariantVector[rownum] <- .5
    } else{
      VariantVector[rownum] <- 0
    }
  }
  Results[,GenomeString] = VariantVector
}


# Variants present in the three low but absent from the high
ConsistentVariants = Results %>% filter(DORN1000==1 &  DORN1088==1 & DORN1038==1 & DORN1194==0 & DORN1037==0 & DORN880==0 & DORN881==0 & DORN976==0 & DORN999==0 & DORN1061==0 & DORN882==0 & DORN933==0 & DORN929==0)

# Variants present in the three low and only partially present in the high 
MediumConsistentVariants= Results %>% filter(DORN1000==1 &  DORN1088==1 & DORN1038==1 & DORN1194<1 & DORN1037<1 & DORN880<1 & DORN881<1 & DORN976<1 & DORN999<1 & DORN1061<1 & DORN882<1 & DORN933<1 & DORN929<1)

# Variants absent from the high but at least partially present in all the low
AlsoMediumConsistentVariants = Results %>% filter(DORN1000>0 &  DORN1088>0 & DORN1038>0 & DORN1194==0 & DORN1037==0 & DORN880==0 & DORN881==0 & DORN976==0 & DORN999==0 & DORN1061==0 & DORN882==0 & DORN933==0 & DORN929==0)


write.csv(Results %>% arrange(Start), file="data/nucdiff925/Comparisons925_1088_All.csv")
