# Amy Campbell
# 2022, updated 2023
# Programmatic comparison of NucDiff results for SA925 vs. each of the other 13 isolates in Patient # 141

library(dplyr)

# Get list of variants & locations in SA1088 compared to SA925
##################################################################

SNPs1088 = read.csv2("Data/WGSData/NucDiffOutput/NucDiff_925_1088/results/SA925_SA1088Comparison_ref_snps.gff",comment.char = "#", sep="\t",header=F)

colnames(SNPs1088) = c("Contig", "NucDiffVersion","SO", "Start", "End", "Blank1", "Blank2", "Blank3", "Description")
SNPs1088 = SNPs1088 %>% select(Start, End, Description)


Structural1088 = read.csv2("Data/WGSData/NucDiffOutput/NucDiff_925_1088/results/SA925_SA1088Comparison_ref_struct.gff",comment.char = "#", sep="\t",header=F)
colnames(Structural1088) = c("Contig", "NucDiffVersion","SO", "Start", "End", "Blank1", "Blank2", "Blank3", "Description")
Structural1088 = Structural1088 %>% select(Start, End, Description)

Additional1088 = read.csv2("Data/WGSData/NucDiffOutput/NucDiff_925_1088/results/SA925_SA1088Comparison_ref_additional.gff",comment.char = "#", sep="\t",header=F)
colnames(Additional1088) = c("Contig", "NucDiffVersion","SO", "Start", "End", "Blank1", "Blank2", "Blank3", "Description")
Additional1088 = Additional1088 %>% select(Start, End, Description)

Variants1088 = rbind(SNPs1088,Structural1088)
Variants1088= rbind(Variants1088, Additional1088)

genomevector=c("SA1088", "SA1000", "SA1038", "SA1194", "SA1061", "SA1037", "SA999", "SA976", "SA933", "SA929", "SA882", "SA881", "SA880")
Results=Variants1088
for(GenomeString in genomevector){
  ComparisonString=paste0("SA925_", GenomeString)
  JustNumbers = stringr::str_remove_all(ComparisonString, "SA")
  print(ComparisonString)
  SNPFile =paste("Data/WGSData/NucDiffOutput/NucDiff_", JustNumbers, "/results/", ComparisonString, "Comparison_ref_snps.gff",sep="")
  AdditionalFile =paste("Data/WGSData/NucDiffOutput/NucDiff_", JustNumbers, "/results/", ComparisonString, "Comparison_ref_additional.gff",sep="")
  StructFile =paste("Data/WGSData/NucDiffOutput/NucDiff_", JustNumbers, "/results/", ComparisonString, "Comparison_ref_struct.gff",sep="")
  
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
ConsistentVariants = Results %>% filter(SA1000==1 &  SA1088==1 & SA1038==1 & SA1194==0 & SA1037==0 & SA880==0 & SA881==0 & SA976==0 & SA999==0 & SA1061==0 & SA882==0 & SA933==0 & SA929==0)

# Variants present in the three low and only partially present in the high 
MediumConsistentVariants= Results %>% filter(SA1000==1 &  SA1088==1 & SA1038==1 & SA1194<1 & SA1037<1 & SA880<1 & SA881<1 & SA976<1 & SA999<1 & SA1061<1 & SA882<1 & SA933<1 & SA929<1)

# Variants absent from the high but at least partially present in all the low
AlsoMediumConsistentVariants = Results %>% filter(SA1000>0 &  SA1088>0 & SA1038>0 & SA1194==0 & SA1037==0 & SA880==0 & SA881==0 & SA976==0 & SA999==0 & SA1061==0 & SA882==0 & SA933==0 & SA929==0)


write.csv(Results %>% arrange(Start), file="Data/WGSData/NucDiffOutput/Comparisons925_1088_All.csv")
