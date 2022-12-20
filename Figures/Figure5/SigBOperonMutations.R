# Amy Campbell
# December 2022 
# Output of aligning each sigB operon gene AA sequence in every genome to that of LAC 
# Where FS and NS mutations were annotated manually w/ tblastN & pgap pseudogene predictions 

library("dplyr")
library("ggplot2")
library("stringr")

# Data loading
###############

CCs = read.csv("data/Phylogeny2022Data/CCMapPlotting.csv")
TreeFilePath = "data/Phylogeny2022Data/RAxML_bestTree.RaxMLTree2022.newick"
STX_data  = read.csv("/Users/amycampbell/Desktop/GriceLabGit/Staphyloxanthin/Data/InVitroData/Phenotypes_Data.csv")
CCs$Genome = CCs$DORN
STX_data$Genome = sapply(STX_data$IsolateID, function(x) str_replace(x, pattern="SA", "DORN") )

TwoPal=RColorBrewer::brewer.pal(8, "Dark2")[c(6,1)]

# RsbU Data
###########
rsbU_mutations= read.csv("/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS/data/RsbU_Mutations.csv")
rsbu_identities = read.csv("/Users/amycampbell/Desktop/GriceLabGit/DFUStrainsWGS/data/rsbU_identities.csv") %>% select(Genome, AA_Identity)
MutationPositionsRsbU = rsbU_mutations %>% left_join(CCs, by="Genome")

# SigB Data
############
sigB_identities= read.csv("data/sigB_Identities.csv") %>% select(Genome, AAIdentity)
sigB_mutations = read.csv("data/SigBMutations.csv")
MismatchPositionsSigB = sigB_mutations %>% left_join(CCs, by="Genome")


# RsbW Data
############
rsbWidentities = read.csv("data/RsbWIdentities.csv") 
rsbW_mutations = read.csv("data/RsbWMismatches.csv")
rsbW_mutations = rsbW_mutations %>% left_join(CCs,by="Genome")


# For 1. RsbU, 2. SigB, 3. RsbW, 4. RsbV:
    # a) Show bar plots of all minor frequency alleles (<110 prevalence) differing non-synonymously (AA change) to LAC 
    # b) Show boxplot for association of having any minor frequency AA change compared to LAC on STX 
    # c) Show comparisons of STX, for subsets of isolates from the same patient belonging to the same CC, of having AA change vs. not (in cases where we can make this comparison)
# Then, show association of having any of these genes pseudogenized with STX production 


##############################
# 1. RSBU ANALYSES  
##############################


#### a) Spatial distribution of of minor frequency, NS mutations compared to LAC
#################################################################################

# get CC assignments from pubMLST 
MutationPositionsRsbU = rsbU_mutations %>% left_join(CCs, by="Genome")

# Check if any are major frequency (>half genomes have it)
print(any(table(MutationPositionsRsbU$NonmatchPos) > 110))

# since not, continue without problems 
MutationBarPlotRsbU = ggplot(MutationPositionsRsbU, aes(x=NonmatchPos, fill=CCLabel)) + geom_bar(width=2) + theme_classic() + scale_fill_brewer(palette="Dark2") + xlab("Mutated AA position (1...333)") + labs(fill="Clonal Complex", y="# Isolates Containing Mutation")
ggsave(MutationBarPlotRsbU, file="data/SigBOperonPlots_STXrevisions/RsbUMutationAA_Barplot.pdf", width=15, height=5)


#### b) Presence/absence of minor frequency, NS mutations compared to LAC vs. STX
##################################################################################

JoinedRsbU = rsbu_identities %>% left_join(STX_data %>%  select(Genome, patient, week_healed, staphyloxanthin), by="Genome") %>% left_join(CCs, by="Genome") 


JoinedRsbU$AnyChanges = if_else(JoinedRsbU$AA_Identity=="Identical", "No", "Yes")
RsbU_mismatches_vs_STX = ggplot(JoinedRsbU, aes(x=AnyChanges, y=log(staphyloxanthin), fill=AnyChanges)) + geom_boxplot() + ggpubr::stat_compare_means(method="t.test") + scale_fill_manual(values=TwoPal) + theme_classic() + xlab("AA Changes in RsbU Compared to LAC's?")
ggsave(RsbU_mismatches_vs_STX, file="data/SigBOperonPlots_STXrevisions/rsbU_identity_vs_STX.pdf",width=5.72, height=4)

#### c) Within-patient, within-CC comparisons of STX in mutant vs. non-mutant isolates for rsbU
################################################################################################
JoinedPositionsRsbU = JoinedRsbU %>% left_join(rsbU_mutations %>% select(Genome, NonmatchPos, ReferenceAA, QueryAA), by="Genome") 

# subset to Patient_CC combinations with 
JoinedPositionsRsbU$patient_CC = paste(JoinedPositionsRsbU$patient, JoinedPositionsRsbU$CCLabel, sep="_")
PatientCCSummary_Rsbu = table(JoinedPositionsRsbU %>% select(AA_Identity, patient_CC))
PatientCC_DF_RsbU = data.frame(PatientCC =names(PatientCCSummary_Rsbu[1,]), NumIdentical = PatientCCSummary_Rsbu[1,], NumNonIdentical= PatientCCSummary_Rsbu[2,])


# Only consider patients who, in at least one CC have at least one 'mutated' isolate and one 'non-mutated'
PatientCCs_With_Combo = PatientCC_DF_RsbU %>% filter(NumIdentical>0 & NumNonIdentical>0)




plotlist_rsbU = list()
i=1
for(patient_combo  in PatientCCs_With_Combo$PatientCC){
  subset_df = (JoinedPositionsRsbU %>% filter(patient_CC==patient_combo))
  if(length(unique(subset_df$NonmatchPos))>2){
    print("Error: More than one mutation in this patient/CC combo!")
  } else{
    PositionMutation = ((subset_df %>% filter(!is.na(NonmatchPos)))$NonmatchPos)[1]
    PatientID = (subset_df$patient)[1]
    
    plot_object = ggplot(subset_df, aes(x=AnyChanges, y=staphyloxanthin, color=AnyChanges)) + geom_jitter(width=.2,height=0, size=3) + ylim(10,100)+
      stat_summary(fun = mean, fun.min = mean, fun.max = mean, geom = "crossbar", width = 0.3, color="black") +
      scale_color_manual(values=TwoPal) + theme_classic() + labs(x="Mutant?", color="Mutant?", y="Staphyloxanthin Production(%502A Absorbance)") + ggtitle( paste0("AA change at AA",toString(PositionMutation), " of RsbU in DFU",toString(PatientID) )) + theme(plot.title=element_text(size=14, hjust=.5, face="bold"))
    plotlist_rsbU[[i]] = plot_object
    i=i+1
    }
  
}

ArrangedRsbUWithinPatient = grid.arrange(plotlist_rsbU[[1]], plotlist_rsbU[[2]], plotlist_rsbU[[3]], plotlist_rsbU[[4]], plotlist_rsbU[[5]], plotlist_rsbU[[6]], ncol=6)
ggsave(ArrangedRsbUWithinPatient, file="data/SigBOperonPlots_STXrevisions/WithinCC_Patient_RsbuXanthin.pdf", width=25,height=5)


##############################
# 2. SIGB ANALYSES  
##############################
# "#00008B" dark blue
# "#06452F" dark green
# "#BF1919" red 
# "#19B7BF" turquoise 
# "#5C0B8F" eggplant 

#### a) Spatial distribution of of minor frequency, NS mutations compared to LAC
#################################################################################
expandedDark2 = RColorBrewer::brewer.pal(8, "Dark2")
expandedDark2 = append(expandedDark2, c("#00008B","#06452F", "#BF1919", "#19B7BF","#5C0B8F" ))

# Are any of the positions not matching LAC's sigB majority frequency alleles? 
print(any(table(MismatchPositionsSigB$NonmatchPos) > 110))

# Filter these out 
MajorityFrequency = names(which(table(MismatchPositionsSigB$NonmatchPos) > 110))
print(MajorityFrequency)
# the mutations at 11 and 256 are both >110 in prevalence 

MismatchPositionsSigB = MismatchPositionsSigB %>% filter(! (as.character(NonmatchPos) %in% MajorityFrequency))

MutationBarPlotSigB= ggplot(MismatchPositionsSigB, aes(x=NonmatchPos, fill=CCLabel)) + geom_bar() + theme_classic() + scale_fill_manual(values=c("#E7298A","#E6AB02", "#5C0B8F")) + xlab("Mutated AA position in LAC SigB") + labs(fill="Clonal Complex")

ggsave(MutationBarPlotSigB, file="data/SigBOperonPlots_STXrevisions/SigBMutationAA_Barplot.pdf", width=7, height=30)



#### b) Presence/absence of minor frequency, NS mutations compared to LAC vs. STX
##################################################################################

JoinedPositionsCCs = STX_data %>% left_join(CCs, by="Genome")
JoinedPositionsSigB = JoinedPositionsCCs %>% left_join(sigB_mutations %>% select(Genome, NonmatchPos, ReferenceAA, QueryAA), by="Genome") 
sigB_identities$Genome = sapply(sigB_identities$Genome, function(x) str_replace(x, pattern="DORN", "SA"))
JoinedPositionsSigB$Genome=JoinedPositionsSigB$IsolateID
JoinedPositionsSigB = JoinedPositionsSigB %>% left_join(sigB_identities, by="Genome")

JoinedPositionsSigB$AnyChanges = if_else(JoinedPositionsSigB$AAIdentity=="Identical", "No", "Yes")

JoinedPositionsSigB = JoinedPositionsSigB %>% filter(!is.na(AAIdentity))


SigB_mismatches_vs_STX = ggplot(JoinedPositionsSigB, aes(x=AnyChanges, y=log(staphyloxanthin), fill=AnyChanges)) + geom_boxplot() + ggpubr::stat_compare_means(method="t.test") + scale_fill_manual(values=TwoPal) + theme_classic() + xlab("AA Changes in SigB Compared to LAC's?")
ggsave(SigB_mismatches_vs_STX, file="data/SigBOperonPlots_STXrevisions/sigB_identity_vs_STX.pdf")


#### c) Within-patient, within-CC comparisons of STX in mutant vs. non-mutant isolates for sigB
################################################################################################

JoinedPositionsSigB$patient_CC = paste(JoinedPositionsSigB$patient, JoinedPositionsSigB$CCLabel, sep="_")
PatientCCSummary_SigB = table(JoinedPositionsSigB %>% select(AAIdentity, patient_CC))


PatientCC_DF_SigB = data.frame(PatientCC =names(PatientCCSummary_SigB[1,]), NumIdentical = PatientCCSummary_SigB[1,], NumNonIdentical= PatientCCSummary_SigB[2,])

# Only consider patients who, in at least one CC have at least one 'mutated' isolate and one 'non-mutated'
# only one in this case! 
PatientCCs_With_Combo = PatientCC_DF_SigB %>% filter(NumIdentical>0 & NumNonIdentical>0)

plotlist_sigB = list()
i=1
for(patient_combo  in PatientCCs_With_Combo$PatientCC){
  subset_df = (JoinedPositionsSigB %>% filter(patient_CC==patient_combo))
  if(length(unique(subset_df$NonmatchPos))>2){
    print("Error: More than one mutation in this patient/CC combo!")
    print(patient_combo)
  } else{
    PositionMutation = ((subset_df %>% filter(!is.na(NonmatchPos)))$NonmatchPos)[1]
    PatientID = (subset_df$patient)[1]
    
    plot_object = ggplot(subset_df, aes(x=AnyChanges, y=staphyloxanthin, color=AnyChanges)) + geom_jitter(width=.2,height=0, size=3) +
      stat_summary(fun = mean, fun.min = mean, fun.max = mean, geom = "crossbar", width = 0.3, color="black") + ylim(0,100)+
      scale_color_manual(values=TwoPal) + theme_classic() + labs(x="Mutant?", color="Mutant?", y="Staphyloxanthin Production(%502A Absorbance)") + ggtitle( paste0("AA change at AA",toString(PositionMutation), " of SigB in DFU",toString(PatientID) )) + theme(plot.title=element_text(size=14, hjust=.5, face="bold"))
    plotlist_sigB[[i]] = plot_object
    i=i+1
  }
  
}


ggsave(plotlist_sigB[[1]], file="data/SigBOperonPlots_STXrevisions/WithinCC_Patient_SigBXanthin.pdf", width=5, height=7)

ArrangedSigBWithinPatient = grid.arrange(plotlist_rsbU[[1]], plotlist_rsbU[[2]], plotlist_rsbU[[3]], plotlist_rsbU[[4]], plotlist_rsbU[[5]], plotlist_rsbU[[6]], ncol=3)
ggsave(ArrangedSigBWithinPatient, file="data/SigBOperonPlots_STXrevisions/WithinCC_Patient_SigBXanthin.pdf", width=20,height=12)



####################
# 2. RSBW ANALYSES  
####################

# Any majority frequency alleles? 
print(any(table(rsbW_mutations$NonmatchPos) > 110))


MutationBarPlotRsbW= ggplot(rsbW_mutations, aes(x=NonmatchPos, fill=CCLabel)) + geom_bar() + theme_classic() + scale_fill_manual(values=c("#1B9E77", "#7570B3","#E7298A" , "#66A61E", "#E6AB02", "#666666", "#5C0B8F")) + xlab("Mutated AA position in LAC RsbW(1...159)") + labs(fill="Clonal Complex")

ggsave(MutationBarPlotRsbW, file="data/SigBOperonPlots_STXrevisions/RsbWMutationAA_Barplot.pdf", width=7, height=5.5)

RsbW_stx = rsbW_mutations %>% left_join(STX_data,by="Genome" )



#### b) Presence/absence of minor frequency, NS mutations compared to LAC vs. STX
##################################################################################

JoinedPositionsCCs = STX_data %>% left_join(CCs, by="Genome")
JoinedPositionsRsbW = JoinedPositionsCCs %>% left_join(rsbW_mutations %>% select(Genome, NonmatchPos, ReferenceAA, QueryAA), by="Genome") 
JoinedPositionsRsbW = JoinedPositionsRsbW %>% left_join(rsbWidentities, by="Genome")

JoinedPositionsRsbW$Genome = sapply(JoinedPositionsRsbW$Genome, function(x) str_replace(x, pattern="DORN", "SA"))


JoinedPositionsRsbW$AnyChanges = if_else(JoinedPositionsRsbW$AAIdentity=="Identical", "No", "Yes")

JoinedPositionsRsbW = JoinedPositionsRsbW %>% filter(!is.na(AAIdentity))


RsbW_mismatches_vs_STX = ggplot(JoinedPositionsRsbW, aes(x=AnyChanges, y=log(staphyloxanthin), fill=AnyChanges)) + geom_boxplot() + ggpubr::stat_compare_means(method="t.test") + scale_fill_manual(values=TwoPal) + theme_classic() + xlab("AA Changes in RsbW Compared to LAC's?")
ggsave(RsbW_mismatches_vs_STX, file="data/SigBOperonPlots_STXrevisions/rsbW_identity_vs_STX.pdf")

#### c) Within-patient, within-CC comparisons of STX in mutant vs. non-mutant isolates for rsbW
################################################################################################

JoinedPositionsRsbW$patient_CC = paste(JoinedPositionsRsbW$patient, JoinedPositionsRsbW$CCLabel, sep="_")
PatientCCSummary_RsbW = table(JoinedPositionsRsbW %>% select(AAIdentity, patient_CC))


PatientCC_DF_RsbW= data.frame(PatientCC =names(PatientCCSummary_RsbW[1,]), NumIdentical = PatientCCSummary_RsbW[1,], NumNonIdentical= PatientCCSummary_RsbW[2,])



# Only consider patients who, in at least one CC have at least one 'mutated' isolate and one 'non-mutated'
# only one in this case! 
PatientCCs_With_Combo = PatientCC_DF_RsbW %>% filter(NumIdentical>0 & NumNonIdentical>0)


plotlist_rsbW = list()
i=1
for(patient_combo  in PatientCCs_With_Combo$PatientCC){
  subset_df = (JoinedPositionsRsbW %>% filter(patient_CC==patient_combo))
  if(length(unique(subset_df$NonmatchPos))>2){
    print("Error: More than one mutation in this patient/CC combo!")
    print(patient_combo)
  } else{
    PositionMutation = ((subset_df %>% filter(!is.na(NonmatchPos)))$NonmatchPos)[1]
    PatientID = (subset_df$patient)[1]
    
    plot_object = ggplot(subset_df, aes(x=AnyChanges, y=staphyloxanthin, color=AnyChanges)) + geom_jitter(width=.2,height=0, size=3) +
      stat_summary(fun = mean, fun.min = mean, fun.max = mean, geom = "crossbar", width = 0.3, color="black") +
      scale_color_manual(values=TwoPal) + theme_classic() + labs(x="Mutant?", color="Mutant?", y="Staphyloxanthin Production(%502A Absorbance)") + ggtitle( paste0("AA change at AA",toString(PositionMutation), " of RsbW in DFU",toString(PatientID) )) + theme(plot.title=element_text(size=14, hjust=.5, face="bold"))
    plotlist_rsbW[[i]] = plot_object
    i=i+1
  }
  
}
plotlist_rsbW[[1]]

plotlist_rsbW[[2]]

ArrangedRsbW = grid.arrange(plotlist_rsbW[[1]],plotlist_rsbW[[2]], ncol=1)
ggsave(ArrangedRsbW, file="data/SigBOperonPlots_STXrevisions/WithinCC_Patient_RsbWXanthin.pdf", height=8, width=5)




ArrangedSigBWithinPatient = grid.arrange(plotlist_rsbU[[1]], plotlist_rsbU[[2]], plotlist_rsbU[[3]], plotlist_rsbU[[4]], plotlist_rsbU[[5]], plotlist_rsbU[[6]], ncol=3)
ggsave(ArrangedSigBWithinPatient, file="data/SigBOperonPlots_STXrevisions/WithinCC_Patient_RsbuXanthin.pdf", width=20,height=12)




















