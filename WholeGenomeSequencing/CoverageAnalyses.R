# Amy Campbell
# Staphyloxanthin paper 2022
# adapted from script to do this to larger # of genomes 
# Make a plot of the coverages by genome, select contigs with coverage outside IQR for blast followup
# Filter to contigs >100bp in length >10X mean coverage

library(dplyr)
library(ggplot2)
library(tidyr)
library(forcats)

coveragestatspath=normalizePath("../Data/WGSData/AllCoverageStats.tsv")
outputfatespath=normalizePath("../Data/WGSData/ContigFates.tsv")

DFU_stats = data.frame(read.table(coveragestatspath, header=T)) 

DFU_statsDORNS = DFU_stats %>% separate(identifer, c("Genome", "ContigID"), "_", remove=F)
DFU_statsDORNS = DFU_statsDORNS %>% group_by(Genome) %>% mutate(MedianContigLength=median(length))


CovgPlot = ggplot(DFU_statsDORNS, aes(x=fct_reorder(as.factor(Genome), mean_cov, median), y=log10(mean_cov+1),  fill=MedianContigLength*.001)) +
  scale_fill_gradient(low="white", high="darkred") + geom_boxplot(outlier.size=.1, outlier.color=NULL, outlier.fill=NULL) +
  xlab("Genome") + ylab("Log10-transformed Mean Coverage by Contig")  + labs(fill="Median Contig Length of Assembly (kb)") +
  theme_bw() + ggtitle("DORN Genome Contig Coverages") +
  theme(legend.position="bottom", axis.text.x=element_text(angle=90), panel.grid.major = element_blank(),
                                                                                                                                                                                                                                                                                                                                                                                                                                                                    panel.grid.minor = element_blank(), legend.title=element_text(size=16), axis.title.x=element_text(size=16), axis.title.y=element_text(size=16), plot.title=element_text(size=20, hjust=.5, face="bold"))  +  geom_hline(yintercept = log10(101), color="#008B8B", linetype="dashed") + geom_hline(yintercept=log10(51),color="#008B8B",linetype="dashed" )
ggsave(CovgPlot, file="CoveragePlot.pdf")

DFU_statsDORNS = DFU_statsDORNS %>% filter(ContigID !="*")

# Filter out contigs <100 bp in length 
print("Prefilter:")
print(dim(DFU_statsDORNS))

DFU_statsDORNS = DFU_statsDORNS %>% filter(length >= 100)

print("Post-length-filter:")
print(dim(DFU_statsDORNS))

# Identify contigs with outlier depths for the genome as falling outside of IQR
DFU_statsDORNS = DFU_statsDORNS %>% group_by(Genome) %>% mutate(LowCutoff = (median(mean_cov)-(1.5*IQR(mean_cov)))) 
DFU_statsDORNS = DFU_statsDORNS %>% group_by(Genome) %>% mutate(HighCutoff = (median(mean_cov)+(1.5*IQR(mean_cov))))

# Filter out contigs with <10X coverage
DFU_statsDORNS = DFU_statsDORNS %>% filter(mean_cov >10)

DFU_statsDORNS = DFU_statsDORNS %>% mutate(isLowOutlier = ((mean_cov < LowCutoff)))
DFU_statsDORNS = DFU_statsDORNS %>% mutate(isHighOutlier = (mean_cov > HighCutoff))

# Mark for followup if outlier in coverage
DFU_statsDORNS = DFU_statsDORNS %>% mutate(SortContig = case_when(isLowOutlier ~ "FollowUp",
                                                                  isHighOutlier ~ "FollowUp",
                                                                  TRUE~"keep"
                                                                  
))

write.csv(DFU_statsDORNS[c("identifer", "Genome", "ContigID", "SortContig")], file=outputfatespath)






