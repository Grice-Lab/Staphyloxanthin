# Amy Campbell
# QC Step: Testing for admixtures
# Analyzing output of the BCF variant calling on reads against cleaned genomes 

library(dplyr)

BCFOutput = read.csv("test_BCF_array.csv")

# Filtration criteria as of 02/12/2021 
# Based in part on https://doi.org/10.1099/mgen.0.000354 
# (Added the overall position depth cutoff since this is the cutoff I've used for my contigs)
# (1) Filter to 10X depth
# (2) Regardless of variant frequency, count SNP sites that are at least 50bp apart from one another
#     (aren't in <50bp 'hot spots' )
# (3) Remove isolates with >30 sites >50bp apart from one another

# this actually doesn't make a difference in the final result in our case
BCFOutput10X = BCFOutput %>% filter(Depth >= 10)

# This function is based directly on the dist_filter() function from
# https://github.com/kumarnaren/mecA-HetSites-calculator/blob/master/filtervcf_v4.py 
# which is what they used in https://doi.org/10.1099/mgen.0.000354
# This is important because this is the metric of '50 bp apart' that they used 
# when benchmarking contamination thresholds 
Filter_Distance50_SNPs = function(positionlist){
  positionlist = sort(positionlist)
  
  j = 1

  # list of unclustered positions
  unclustered = c()
  
  # list of clusters
  clusters = c()
  
  currentcluster = c()
  if(length(positionlist) == 1){
    return(positionlist)
  }else{

    for(i in 1:(length(positionlist) - 1)){
      if(positionlist[i+1] - positionlist[i] < 50){
        currentcluster = append(currentcluster, c(positionlist[i], positionlist[i+1]))
        clusters[[j]] = currentcluster
      } else{
        unclustered = append(unclustered,c(positionlist[i], positionlist[i+1]))
        j = j + 1
        currentcluster = c()
      }
    }
    totalClustered = unlist(clusters)
    filteredSNPs = setdiff(unclustered,totalClustered)
    return(filteredSNPs) }
}

# Count >50 bp apart SNPs 
k=0
over30array =c()
for(genome in unique(BCFOutput10X$X..Genome)){
  positionlist = (BCFOutput10X %>% filter(X..Genome == genome))$Position
  SNPlist = Filter_Distance50_SNPs(positionlist)
  print(length(SNPlist))
  if(length(SNPlist) > 30 ){
    print(length(SNPlist))
    over30array = append(over30array, genome)
    k=k+1
  }
}
print(over30array)

DFUisolatesInfo = read.csv("DFU_Staph_aureus_isolates.csv")
DFUisolatesInfo$DORN = paste0("DORN", DFUisolatesInfo$Doern.lab.bank.)

DFUisolatesInfo = DFUisolatesInfo %>% select(DORN, patient_id,visit)
DFUisolatesInfo_Subset = DFUisolatesInfo %>% filter(DORN %in% over30array)
write.csv(DFUisolatesInfo_Subset, file = "ContaminatedIsolates02-03-22.csv")

