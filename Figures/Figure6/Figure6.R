library("ggplot2")
library("gridExtra")
library("dplyr")
library("gridExtra")

# Set this to wherever you keep the repo. 'Data' must be 
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


# Mean of each normalized phenotype by 
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

