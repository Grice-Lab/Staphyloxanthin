# Amy Campbell
# 07/2022 
# Modified from larger 04/2022 script (DEanalysis.R) to specifically reproduce figures included in Figure 5
# Differential expression analysis for SA925 and SA1088 exposed to H202 in vitro 

library(dplyr)
library(Rsubread)
library(Rsamtools)
library(stringr)
library(DESeq2)
library(ggplot2)
library(tidyr)
library(viridis)

setwd("/Users/amycampbell/Desktop/GriceLabGit/Staphyloxanthin/")


SampleMapPath = "Data/RNASeq/MappingRNASeq.txt"
GeneNameMapPath = "Data/RNASeq/GeneNamesPreliminary.csv"
RNACountsPath = "Data/RNASeq/RNACounts.rda"
bamdirpath = "Data/RNASeq/bams/"
GTFPath = "Data/RNASeq/SA925_ref.gff"
StressGeneMap = "Data/RNASeq/StressResponseGenesInclude.csv"
PhageCDS=read.csv("Data/RNASeq/PhageCDS.csv") %>% select(PGAPid)


########################
# 1.) Get FeatureCounts
########################

bamdirectory = list.files(bamdirpath)
bamdirectory = sapply(bamdirectory, function(x) paste0(bamdirpath, x))
bamfiles=BamFileList(bamdirectory)

gtffile <- file.path(GTFPath)

## Uncomment if you want to specifically reproduce the featurecounts step (takes some time is why it's commented out; otherwise, just load in next step)
# MyCounts = featureCounts(files = bamdirectory,annot.ext =gtffile, isGTFAnnotationFile=T, GTF.featureType="CDS", GTF.attrType = "ID", isPairedEnd=T)
# save(MyCounts, file="/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS/RNAseqInputs/RNACounts.rda")


#################################
# 2.) Load in, process other data
#################################

# To skip the above step, load this so you get the MyCounts variable
load(RNACountsPath) 

samplemap=read.csv2(SampleMapPath, sep="\t", header=F)

# Read in sample map and gene map 
samplemap=read.csv2(SampleMapPath, sep="\t", header=F)
colnames(samplemap) = c("Sample", "SampleTreatment", "Condition")
genemap=read.csv2(GeneNameMapPath, sep=',')
colnames(samplemap) = c("Sample", "SampleTreatment", "Condition")


#  Format genemap
################# 
genemap = genemap %>% mutate(GeneAdjust= if_else( (is.na(Gene) | Gene==""), RefSeqID, Gene))

# fixing a few annotations so their identities are clearer 
genemap$GeneAdjust = if_else(genemap$GeneAdjust=="SAOUHSC_01999", "bcp", genemap$GeneAdjust)
genemap$GeneAdjust = if_else(genemap$AnnotID=="cds-pgaptmp_002613", "mgrA", genemap$GeneAdjust)
genemap$GeneAdjust = if_else(genemap$AnnotID=="cds-pgaptmp_001532", "mcsA", genemap$GeneAdjust)
genemap$GeneAdjust = if_else(genemap$GeneAdjust=="SAOUHSC_02171", "sak", genemap$GeneAdjust)



# Make a plot of Sample vs. Total Reads assigned; 1 replicate
# of SA925_control has the highest # reads by an order of magnitude; 
# 1 replicate of SA925_Treatment failed during sequencing so was excluded
#########################################################################
CountsBySample=data.frame((colSums(MyCounts$counts)))
CountsBySample$Sample=sapply(row.names(CountsBySample), function(x) str_split(x, pattern="_")[[1]][1])
colnames(CountsBySample) = c("TotalAssignedReads", "Sample")

CountsBySample = samplemap %>% left_join(CountsBySample, by="Sample")
CountsBySample[is.na(CountsBySample)] <- 0
CountsBySample$Strain = sapply(CountsBySample$SampleTreatment, function(x) if_else(grepl(x, pattern="ARM72"),"SA925(high)","SA1088(low)") )

CountsBySample$ConditionStrain = paste(CountsBySample$Strain, CountsBySample$Condition, sep=" : ")
CountsBySample$ConditionStrain = factor(CountsBySample$ConditionStrain )
TotalCountsPlot = ggplot(CountsBySample, aes(x=Sample, y=TotalAssignedReads, fill=ConditionStrain)) + geom_bar(stat="identity", color="black")
TotalCountsPlot$data$Sample = factor(TotalCountsPlot$data$Sample , levels=(TotalCountsPlot$data %>% arrange(ConditionStrain))$Sample)

TotalCountsPlot$data$Sample = factor(TotalCountsPlot$data$Sample , levels=(TotalCountsPlot$data %>% arrange(Strain, Condition))$Sample)

treatmentcombos = c("#BCD2E8", "#2E5984", "#FCF787", "#EAAA00")
TotalCountsPlot=TotalCountsPlot + scale_fill_manual(values=treatmentcombos) + scale_y_continuous(limits=c(0, 25000000), breaks=seq(0, 25000000, 1000000))


# Read counts data into DEseq
##############################
SampleMapping = CountsBySample %>% filter(Sample!="ARM6")
CountsBySample$TotalAssignedReads = NULL

columndataDE = SampleMapping %>% select(Condition, Strain)
columndataDE$Condition= factor(columndataDE$Condition)

columndataDE$Strain = sapply(columndataDE$Strain, function(x) str_split(x, pattern='\\(')[[1]][1])
columndataDE$Strain= factor(columndataDE$Strain)

rownames(columndataDE) = SampleMapping$Sample

CountsMatrix = MyCounts$counts
colnames(CountsMatrix) = sapply(colnames(CountsMatrix), function(x) str_split(x, pattern="_")[[1]][1])
CountsMatrix = CountsMatrix[,row.names(columndataDE)]

# "group" is the condition_strain (so ctrl_SA925, treatment_925, etc.)
columndataDE$group = factor(paste(columndataDE$Condition, columndataDE$Strain, sep="_"))

DEData =DESeqDataSetFromMatrix(countData=CountsMatrix, 
                               colData=columndataDE,
                               design= ~group)
myDESeqObj = DESeq(DEData)


# PCA on VST-transformed
########################
myDESeqObjVST = vst(myDESeqObj)
plotPCA(myDESeqObjVST,intgroup=c("group")) + scale_color_manual(values=c("#BCD2E8", "#2E5984", "#FCF787", "#EAAA00"))

# Checking mapping:sak should be 0 counts in SA1088 samples since there's no sak gene in 1088
##########################################################################################
CountsMatrix[c("cds-pgaptmp_000025"),]
# 0 in all of what are supposed to be 1088 genomes so that's good. (ARMs 3,4, 7,8, 11, 12)


#############################################################################
# 3.) Prepare dataframes for the following comparisons:
#     (i)  results_SA925_StressDF -- SA925 with H202 vs. SA925 without 
#     (ii) results_SA1088_StressDF -- SA1088 with H202 vs. SA1088 without 
#     (iii) results_genotypeDF -- SA925 vs. SA1088 under control conditions
#     (iv) results_genotypeH202DF -- SA925 vs. SA1088 under H202 conditions
#############################################################################

# (i)
results_SA925_Stress <- results(myDESeqObj, contrast=c("group", "treatment_SA925","ctrl_SA925"))
results_SA925_StressDF = (data.frame(results_SA925_Stress))
results_SA925_StressDF$AnnotID = row.names(results_SA925_StressDF)
results_SA925_StressDF = results_SA925_StressDF %>% left_join(genemap, by="AnnotID")

# (ii)
results_SA1088_Stress <- results(myDESeqObj, contrast=c("group", "treatment_SA1088","ctrl_SA1088"))
results_SA1088_StressDF = (data.frame(results_SA1088_Stress))
results_SA1088_StressDF$AnnotID = row.names(results_SA1088_StressDF)
results_SA1088_StressDF = results_SA1088_StressDF %>% left_join(genemap, by="AnnotID")

# (iii)
results_genotype <- results(myDESeqObj, contrast=c("group", "ctrl_SA925","ctrl_SA1088"))
results_genotypeDF = data.frame(results_genotype) %>% select(padj,pvalue, log2FoldChange)
results_genotypeDF$AnnotID = row.names(results_genotypeDF)
results_genotypeDF = results_genotypeDF %>% left_join(genemap, by="AnnotID")


# (iv)
results_genotypeH202 <- results(myDESeqObj, contrast=c("group", "treatment_SA925","treatment_SA1088"))
results_genotypeH202DF = data.frame(results_genotypeH202) %>% select(padj,pvalue, log2FoldChange)
results_genotypeH202DF$AnnotID = row.names(results_genotypeH202DF)
results_genotypeH202DF = results_genotypeH202DF %>% left_join(genemap, by="AnnotID")



######################################################################
# 4. Heatmaps of common stress response genes, selected DE toxin genes
######################################################################

stressGenes = read.csv(StressGeneMap)

# Selected subset of toxin-related genes which are DE between strains 
AnnotIDsToxins = c("cds-pgaptmp_000020", "cds-pgaptmp_000025", "cds-pgaptmp_000027",
                   "cds-pgaptmp_001239", "cds-pgaptmp_001240", "cds-pgaptmp_001993",
                   "cds-pgaptmp_002323", "cds-pgaptmp_002324", "cds-pgaptmp_002725",
                   "cds-pgaptmp_002743", "cds-pgaptmp_002745", "cds-pgaptmp_001629", 
                   "cds-pgaptmp_001630", "cds-pgaptmp_001631","cds-pgaptmp_001632",
                   "cds-pgaptmp_000173", "cds-pgaptmp_000939")

# Stress genes
###############

results_SA925_StressDFjoined = results_SA925_StressDF %>% left_join(stressGenes, by="AnnotID")
StressResponsive925stress = results_SA925_StressDFjoined %>% filter(!is.na(GeneName))
StressResponsive925stress = StressResponsive925stress %>% group_by(Grouping, OperonRegulonSubgrouping) %>% arrange(Grouping,OperonRegulonSubgrouping,  GeneName) %>% ungroup()

StressResponsive925stress = StressResponsive925stress %>% mutate(GeneAdjust= if_else(is.na(GeneName), GeneAdjust, GeneName))
StressResponsive925stress$comparison = "SA925_h2o2"
StressResponsive925 = StressResponsive925stress %>% select(GeneAdjust, comparison, log2FoldChange, padj)
order = StressResponsive925stress$GeneAdjust

# Stress genes in SA1088 re: stress
results_SA1088_Stressjoined = results_SA1088_StressDF %>% left_join(stressGenes, by="AnnotID")
results_SA1088stress = results_SA1088_Stressjoined %>% filter(!is.na(GeneName))
results_SA1088stress = results_SA1088stress %>% group_by(Grouping, OperonRegulonSubgrouping) %>% arrange(Grouping,OperonRegulonSubgrouping,  GeneName) %>% ungroup()
StressResponsive1088stress = results_SA1088stress %>% mutate(GeneAdjust= if_else(is.na(GeneName), GeneAdjust, GeneName))
StressResponsive1088stress$comparison = "SA1088_h2o2"
StressResponsive1088stress = StressResponsive1088stress %>% select(GeneAdjust, comparison, log2FoldChange, padj)

# Stress genes in 925 vs. 1088 under H202 conditions
results_Genotype_joinedH202 = results_genotypeH202DF %>% left_join(stressGenes, by="AnnotID")
results_Genotype_joinedH202 = results_Genotype_joinedH202 %>% filter(!is.na(GeneName))
results_Genotype_joinedH202 = results_Genotype_joinedH202 %>% group_by(Grouping, OperonRegulonSubgrouping) %>% arrange(Grouping,OperonRegulonSubgrouping,  GeneName) %>% ungroup()
StressResponsiveGenotypeH202 = results_Genotype_joinedH202 %>% mutate(GeneAdjust= if_else(is.na(GeneName), GeneAdjust, GeneName))
StressResponsiveGenotypeH202$comparison = "SA925_1088_H202"
StressResponsiveGenotypeH202 = StressResponsiveGenotypeH202 %>% select(GeneAdjust, comparison, log2FoldChange, padj)



StressGenesAll = rbind(StressResponsive925, StressResponsive1088stress)
StressGenesAll = rbind(StressGenesAll, StressResponsiveGenotypeH202)

StressGenesMelted = StressGenesAll %>% reshape2::melt(id.vars=c("comparison", "GeneAdjust"))
StressGenesMeltedPvalues =  StressGenesMelted %>% filter(variable=="padj")

StressGenesMeltedLFCs =  StressGenesMelted %>% filter(variable=="log2FoldChange")

StressGenesMeltedPvalues = StressGenesMeltedPvalues %>%  mutate(sigLabel=case_when(value < .05 & value >= .01 ~"*",
                                                                                   value <.01 & value >= .001 ~"**",
                                                                                   value < .001 & value >= .0001 ~ "***",
                                                                                   value < .0001~ "****",
                                                                                   
                                                                                   TRUE~ ""))
StressGenesMeltedPvalues
StressGenesMeltedLFCs = StressGenesMeltedLFCs %>% left_join(StressGenesMeltedPvalues %>% select(comparison, GeneAdjust, sigLabel), by=c("comparison", "GeneAdjust"))
StressGenesMeltedLFCs


StressPlot = ggplot(StressGenesMeltedLFCs, aes(x=comparison, y=GeneAdjust, fill=value)) + geom_tile() + scale_fill_viridis(option="plasma") + geom_text(aes(label=sigLabel), vjust=.7, color="white") + theme_classic() + labs(fill="log2FC") + ylab("Gene")
StressPlot$data$GeneAdjust = factor(StressPlot$data$GeneAdjust, levels=rev(order))
StressPlot$data$comparison = factor(StressPlot$data$comparison, levels=c("SA925_1088_H202", "SA925_h2o2", "SA1088_h2o2"))

grepl(Pattern, lst_B)
ggsave(StressPlot, file="/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS/XanthinPaperFigures/UTDheatmap_stressOnly.pdf", height=15, width=8)

# Filter to significant genes, or genes in operons with significant genes
#########################################################################
# operon strings to keep full operon of
signif_operons_keep_strings=c("wal","mut","pol", "rec", "ssb", "add","uvr","czr", "mnt", "sbn", "sir","sfa","clp","mcs","suf","gro")
patterntest=paste(signif_operons_keep_strings,collapse="|")

# Phage CDS
###########
PhageCDS$AnnotID = sapply(PhageCDS$PGAPid, function(x) str_remove(x,"ID="))

StressGenesInPhage = stressGenes %>% filter(AnnotID %in% PhageCDS$AnnotID )

StressGenesMeltedLFCsKeepGenes= StressGenesMeltedLFCs %>% filter( (sigLabel %in% c("*", "**", "***", "****")) | grepl(patterntest, GeneAdjust))
StressGenesMeltedLFCs_Filtered = StressGenesMeltedLFCs %>% filter(GeneAdjust %in% StressGenesMeltedLFCsKeepGenes$GeneAdjust )
StressGenesMeltedLFCs_Filtered = StressGenesMeltedLFCs_Filtered %>% filter(!(GeneAdjust %in% StressGenesInPhage$GeneName))
StressPlotFiltered = ggplot(StressGenesMeltedLFCs_Filtered, aes(x=comparison, y=GeneAdjust, fill=value)) + geom_tile() + scale_fill_viridis(option="plasma") + geom_text(aes(label=sigLabel), vjust=.7, color="white") + theme_classic() + labs(fill="log2FC") + ylab("Gene")
StressPlotFiltered$data$GeneAdjust = factor(StressPlotFiltered$data$GeneAdjust, levels=rev(order))
StressPlotFiltered$data$comparison = factor(StressPlotFiltered$data$comparison, levels=c("SA925_1088_H202", "SA925_h2o2", "SA1088_h2o2"))

ggsave(StressPlotFiltered, file="Figures/Figure5/UTD_Heatmap.pdf", height=15, width=5.75)

# Making heatmap with select toxins
######################################
AnnotIDsToxins = c("cds-pgaptmp_000020", "cds-pgaptmp_000025", "cds-pgaptmp_000027",
                   "cds-pgaptmp_001239", "cds-pgaptmp_001240", "cds-pgaptmp_001993",
                   "cds-pgaptmp_002323", "cds-pgaptmp_002324", "cds-pgaptmp_002725",
                   "cds-pgaptmp_002743", "cds-pgaptmp_002745", "cds-pgaptmp_001629", 
                   "cds-pgaptmp_001630", "cds-pgaptmp_001631","cds-pgaptmp_001632", "cds-pgaptmp_000173", "cds-pgaptmp_000939")
# pgaptmp_000020 -- entA
# pgaptmp_000025 -- sak
# pgaptmp_000027 -- scn
# pgaptmp_001239 -- seL/entL 
# pgaptmp_001240 -- sec2/entC2
# pgaptmp_001993 -- entH
# pgaptmp_002323 -- hlgB
# pgaptmp_002324 -- hlgC
# pgaptmp_002725 -- hld
# pgaptmp_002743 -- hlb
# pgaptmp_002745 -- entK (not EntQ)
# pgaptmp_001629 PSmAlpha-1
# pgaptmp_001630 PSMalpha-2
# pgaptmp_001631 PSMalpha-3
# pgaptmp_001632 PSMalpha-4
# pgaptmp_000173 -- traP
# pgaptmp_000939 -- hly





Toxin925 = results_SA925_StressDF %>% filter(AnnotID %in% AnnotIDsToxins) %>% select(GeneAdjust, log2FoldChange, padj) 
Toxin925$Comparison = "SA925_Stress"

Toxin1088 = results_SA1088_StressDF %>% filter(AnnotID %in% AnnotIDsToxins) %>% select(GeneAdjust, log2FoldChange, padj) 
Toxin1088$Comparison = "SA1088_Stress" 

ToxinGenotypeControl = results_genotypeDF %>% filter(AnnotID %in% AnnotIDsToxins) %>% select(GeneAdjust, log2FoldChange, padj) 
ToxinGenotypeControl$Comparison = "SA925_SA1088_Control"

ToxinGenotypeStress = results_genotypeH202DF %>% filter(AnnotID %in% AnnotIDsToxins) %>% select(GeneAdjust, log2FoldChange, padj) 
ToxinGenotypeStress$Comparison = "SA925_SA1088_Stress"


MeltedToxins =do.call(rbind, list(Toxin925, Toxin1088,ToxinGenotypeControl, ToxinGenotypeStress ))
MeltedToxins$GeneAdjust = if_else(MeltedToxins$GeneAdjust=="SAOUHSC_02171", "sak", MeltedToxins$GeneAdjust)

MeltedToxins = MeltedToxins %>% mutate(ptext = case_when(padj < .05 & padj >= .01 ~"*", 
                                                         padj <.01 & padj >= .001 ~"**",
                                                         padj < .001 & padj >= .0001 ~ "***", 
                                                         padj < .0001~ "****",
                                                         
                                                         TRUE~ ""))


ToxinHeatMap = ggplot(MeltedToxins, aes(x=Comparison, y=GeneAdjust, fill=log2FoldChange)) + geom_tile() + scale_fill_viridis(option="plasma")  + geom_text(aes(label=ptext),color="white")+ theme_classic() +theme(axis.text.y=element_text(size=15)) + labs(fill="log2FC") + ylab("Gene")
ToxinHeatMap$data$Comparison = factor(ToxinHeatMap$data$Comparison, levels=c("SA925_SA1088_Control", "SA925_SA1088_Stress", "SA925_Stress", "SA1088_Stress"))



