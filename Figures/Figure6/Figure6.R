library("ggplot2")
library("gridExtra")
library("dplyr")
library("gridExtra")
library("bestNormalize")
library("stringr")

# Set this to wherever you keep the repo. 'Data' must be in this folder 
setwd("/Users/amycampbell/Desktop/GriceLabGit/Staphyloxanthin/")

FullData = read.csv("Data/InVitroData/staphyloxanthin_paper_updated.csv")

FullData$Healed = if_else(is.na(FullData$week_healed ), "Unhealed", "Healed")


# Selecting normalization methods for each phenotype 
# (ARM)
####################################################

# Staphyloxanthin
stx <-  FullData$staphyloxanthin
log_stx <-  log(FullData$staphyloxanthin)
shapiro.test(stx)
shapiro.test(log_stx)
FullData$logStaphyloxanthin = log_stx
# log() selected

# biofilm 
biofilm <- FullData$biofilm
shapiro.test(biofilm)
log_biofilm<-log(biofilm)
shapiro.test(log_biofilm)
FullData$logBiofilm = log_biofilm
# log() selected


# siderophore(using bestNormalize package bc log doesnt render normal)
siderophore=FullData$siderophore
shapiro.test(FullData$siderophore)
shapiro.test(log(FullData$siderophore))
bestNormalize(FullData$siderophore)
siderophore_ordernorm = predict(orderNorm(FullData$siderophore))
shapiro.test(siderophore_ordernorm)
FullData$siderophore_normalized = siderophore_ordernorm


# staphylokinase (once again using bestNormalize package)
kinasedata = FullData$staphylokinase
shapiro.test(kinasedata)
shapiro.test(log(kinasedata))
bestNormalize(kinasedata)
staphylokianse_ordernorm = predict(orderNorm(kinasedata))
shapiro.test(staphylokianse_ordernorm)

FullData$staphylokinase_normalized = staphylokianse_ordernorm

#####################################################
# Figure 6B. Correlation between STX in vitro and time
############################3#########################
# Mean of each normalized phenotype by patient/visit combination
MeanByVisit = FullData %>% group_by(patient,visit) %>% summarize(patient=patient, visit=visit,
                                                                 staphyloxanthinmean=mean(logStaphyloxanthin),
                                                                 biofilmmean=mean(logBiofilm),
                                                                 siderophoremean=mean(siderophore_normalized), 
                                                                 staphylokinasemean=mean(staphylokinase_normalized, na.rm=T),
                                                                 week_healed=week_healed) %>% unique()


MeanByVisit_healed = MeanByVisit %>% filter(!is.na(week_healed))
MeanByVisit_unhealed = MeanByVisit %>% filter(is.na(week_healed))

cor.test(MeanByVisit_unhealed$visit, log(MeanByVisit_unhealed$staphyloxanthinmean))
cor.test(MeanByVisit_healed$visit, log(MeanByVisit_healed$staphyloxanthinmean))

HealedPlotTime_STX = ggplot(MeanByVisit_healed, aes(x=visit, y=staphyloxanthinmean)) + geom_point(color="dodgerblue") + geom_smooth(method="lm") + theme_classic() #(lm(MeanByVisit_healed$staphyloxanthinmean~ MeanByVisit_healed$visit))
UnhealedPlotTime_STX = ggplot(MeanByVisit_unhealed, aes(x=visit, y=staphyloxanthinmean)) + geom_point(color="darkorange") + geom_smooth(method="lm") + theme_classic() #(lm(MeanByVisit_healed$staphyloxanthinmean~ MeanByVisit_healed$visit))
ggsave(grid.arrange(HealedPlotTime_STX,UnhealedPlotTime_STX , ncol=2), file="Figures/Figure6/CorrelationTimeSTX.pdf", width=8, height=4)







TwoPal=RColorBrewer::brewer.pal(8, "Dark2")[c(6,1)]

WeekPalette =RColorBrewer::brewer.pal(11, "Spectral")

darkpurple =RColorBrewer::brewer.pal(11, "PuOr")[10]
WeekPalette = c("#8B0000", WeekPalette)

WeekPalette = c("black", WeekPalette)
WeekPalette = append(WeekPalette, darkpurple)
colormap=data.frame(colors=WeekPalette, week=c(0,2,4,6,8,10,12,14,16,18,20,22,24,26))

#################################################################
# Figure 6A. Within-patient, within-CC sigma B operon comparisons
##################################################################
# CC assignments 
CCassignments = read.csv("Data/SigBOperonSequences/CCMapPlotting.csv")

PaperData = FullData
PaperData$Genome=PaperData$IsolateID
MappingPatient_Timepoint = PaperData %>% select(patient, visit,IsolateID)
PaperData$week = factor(2*PaperData$visit)
# RsbU Data
###########
rsbU_mutations= read.csv("Data/SigBOperonSequences/RsbU_Mutations.csv")
rsbu_identities = read.csv("Data/SigBOperonSequences/rsbU_identities.csv") %>% select(Genome, AA_Identity)
MutationPositionsRsbU = rsbU_mutations %>% left_join(CCassignments, by="Genome")

# SigB Data
############
sigB_identities= read.csv("Data/SigBOperonSequences/sigB_Identities.csv") %>% select(Genome, AAIdentity)
sigB_mutations = read.csv("Data/SigBOperonSequences/SigBMutations.csv")
MismatchPositionsSigB = sigB_mutations %>% left_join(CCassignments, by="Genome")

# RsbW Data
############
rsbWidentities = read.csv("Data/SigBOperonSequences/RsbWIdentities.csv") 
rsbW_mutations = read.csv("Data/SigBOperonSequences/RsbWMismatches.csv")
rsbW_mutations = rsbW_mutations %>% left_join(CCassignments,by="Genome")


# RsbU NS mutations 
###################

# Check if any are major frequency (>half genomes have it)
print(any(table(MutationPositionsRsbU$NonmatchPos) > 110))


JoinedRsbU = rsbu_identities %>% left_join(PaperData %>%  select(Genome, patient, week, week_healed, staphyloxanthin), by="Genome") %>% left_join(CCassignments, by="Genome") 
JoinedRsbU$AnyChanges = if_else(JoinedRsbU$AA_Identity=="Identical", "No", "Yes")

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
  cc = str_split(patient_combo, "_")[[1]][2]
  if(length(unique(subset_df$NonmatchPos))>2){
    print("Error: More than one mutation in this patient/CC combo!")
  } else{
    PositionMutation = ((subset_df %>% filter(!is.na(NonmatchPos)))$NonmatchPos)[1]
    PatientID = (subset_df$patient)[1]
    colorpal= colormap %>% filter(week %in% subset_df$week )
    plot_object = ggplot(subset_df, aes(x=AnyChanges, y=staphyloxanthin, color=week)) + geom_jitter(width=.2,height=0, size=3) + ylim(10,100)+
      stat_summary(fun = mean, fun.min = mean, fun.max = mean, geom = "crossbar", width = 0.3, color="black") + 
      scale_color_manual(values=colorpal$colors) + theme_classic() + labs(x="Mutant?", color="Week Collected", y="Staphyloxanthin Production(%502A Absorbance)") + ggtitle( paste0("AA change at AA",toString(PositionMutation), " of RsbU in DFU",toString(PatientID) ,"(",cc,")") ) + theme(plot.title=element_text(size=14, hjust=.5, face="bold"))
    plot_object$data$week = factor(plot_object$data$week , levels=c(0,2,4,6,8,10,12,14,16,18,20,22,24,26))
    plotlist_rsbU[[i]] = plot_object
    i=i+1
  }
  
}

ArrangedRsbUWithinPatient = grid.arrange(plotlist_rsbU[[1]], plotlist_rsbU[[2]], plotlist_rsbU[[3]], plotlist_rsbU[[4]], plotlist_rsbU[[5]], plotlist_rsbU[[6]], ncol=6)
ggsave(ArrangedRsbUWithinPatient, file="Figures/Figure6/WithinCC_Patient_RsbuXanthin.pdf", width=25,height=5)

# SigB 
#################################

JoinedPositionsCCs = PaperData %>% left_join(CCassignments, by="Genome")
JoinedPositionsSigB = JoinedPositionsCCs %>% left_join(sigB_mutations %>% select(Genome, NonmatchPos, ReferenceAA, QueryAA), by="Genome") 
sigB_identities$Genome = sapply(sigB_identities$Genome, function(x) str_replace(x, pattern="DORN", "SA"))
JoinedPositionsSigB$Genome=JoinedPositionsSigB$IsolateID
JoinedPositionsSigB = JoinedPositionsSigB %>% left_join(sigB_identities, by="Genome")

JoinedPositionsSigB$AnyChanges = if_else(JoinedPositionsSigB$AAIdentity=="Identical", "No", "Yes")

JoinedPositionsSigB = JoinedPositionsSigB %>% filter(!is.na(AAIdentity))



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
for(patient_combo in PatientCCs_With_Combo$PatientCC){
  cc = str_split(patient_combo, "_")[[1]][2]
  
  subset_df = (JoinedPositionsSigB %>% filter(patient_CC==patient_combo))
  if(length(unique(subset_df$NonmatchPos))>2){
    print("Error: More than one mutation in this patient/CC combo!")
    print(patient_combo)
  } else{
    PositionMutation = ((subset_df %>% filter(!is.na(NonmatchPos)))$NonmatchPos)[1]
    PatientID = (subset_df$patient)[1]
    colorpal= colormap %>% filter(week %in% subset_df$week )
    plot_object = ggplot(subset_df, aes(x=AnyChanges, y=staphyloxanthin, color=week)) + geom_jitter(width=.2,height=0, size=3) +
      stat_summary(fun = mean, fun.min = mean, fun.max = mean, geom = "crossbar", width = 0.3, color="black") + ylim(0,100)+
      scale_color_manual(values=colorpal$colors) + theme_classic() + labs(x="Mutant?", color="Week collected", y="Staphyloxanthin Production(%502A Absorbance)") + ggtitle( paste0("AA change at AA",toString(PositionMutation), " of SigB in DFU",toString(PatientID), "(", cc,")" )) + theme(plot.title=element_text(size=14, hjust=.5, face="bold"))
    plotlist_sigB[[i]] = plot_object
    i=i+1
  }
  
}


ggsave(plotlist_sigB[[1]], file="Figures/Figure6/WithinCC_Patient_SigBXanthin.pdf", width=4, height=5)

ArrangedSigBWithinPatient = grid.arrange(plotlist_sigB[[1]], ncol=3)
ggsave(ArrangedSigBWithinPatient, file="Figures/Figure6/WithinCC_Patient_SigBXanthin.pdf", width=20,height=12)


# RsbW
###############
# Any majority frequency alleles? 
print(any(table(rsbW_mutations$NonmatchPos) > 110))

RsbW_stx = rsbW_mutations %>% left_join(PaperData,by="Genome" )

JoinedPositionsCCs = PaperData %>% left_join(CCassignments, by="Genome")
JoinedPositionsRsbW = JoinedPositionsCCs %>% left_join(rsbW_mutations %>% select(Genome, NonmatchPos, ReferenceAA, QueryAA), by="Genome") 
JoinedPositionsRsbW = JoinedPositionsRsbW %>% left_join(rsbWidentities, by="Genome")

JoinedPositionsRsbW$Genome = sapply(JoinedPositionsRsbW$Genome, function(x) str_replace(x, pattern="DORN", "SA"))


JoinedPositionsRsbW$AnyChanges = if_else(JoinedPositionsRsbW$AAIdentity=="Identical", "No", "Yes")

JoinedPositionsRsbW = JoinedPositionsRsbW %>% filter(!is.na(AAIdentity))

#Within-patient, within-CC comparisons of STX in mutant vs. non-mutant isolates for rsbW
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
  cc = str_split(patient_combo, "_")[[1]][2]
  
  if(length(unique(subset_df$NonmatchPos))>2){
    print("Error: More than one mutation in this patient/CC combo!")
    print(patient_combo)
  } else{
    PositionMutation = ((subset_df %>% filter(!is.na(NonmatchPos)))$NonmatchPos)[1]
    PatientID = (subset_df$patient)[1]
    colorpal= colormap %>% filter(week %in% subset_df$week )
    
    plot_object = ggplot(subset_df, aes(x=AnyChanges, y=staphyloxanthin, color=week)) + geom_jitter(width=.2,height=0, size=2) +
      stat_summary(fun = mean, fun.min = mean, fun.max = mean, geom = "crossbar", width = 0.3, color="black") + ylim(0,100)+
      scale_color_manual(values=colorpal$colors) + theme_classic() + labs(x="Mutant?", color="Week collected", y="Staphyloxanthin Production(%502A Absorbance)") + ggtitle( paste0("AA change at AA",toString(PositionMutation), " of RsbW in DFU",toString(PatientID) )) + theme(plot.title=element_text(size=14, hjust=.5, face="bold"))
    plotlist_rsbW[[i]] = plot_object
    i=i+1
  }
  
}
plotlist_rsbW[[1]]

plotlist_rsbW[[2]]

ArrangedRsbW = grid.arrange(plotlist_rsbW[[1]],plotlist_rsbW[[2]], ncol=2)
ggsave(ArrangedRsbW, file="Figures/Figure6/WithinCC_Patient_RsbWXanthin.pdf", height=4, width=7)





