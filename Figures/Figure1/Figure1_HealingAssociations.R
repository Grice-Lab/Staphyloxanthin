# Figure 1
# Healing associations
# Modified by Amy Campbell from Amelia McCready-Vangi's and Amy Campbell's scripts
# To reproduce figures 1B-1I

library("ggplot2")
library("gridExtra")
library("lme4")
library("dplyr")
library("bestNormalize")

set.seed(19104)
setwd("/Users/amycampbell/Desktop/GriceLabGit/Staphyloxanthin/")

FullData = read.csv("Data/InVitroData/staphyloxanthin_paper_updated.csv")

FullData$Healed = if_else(is.na(FullData$week_healed ), "Unhealed", "Healed")
FullData$Healed_by_12 = if_else(is.na(FullData$week_healed ) | FullData$week_healed > 12, "Unhealed12", "Healed12")

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

# Figure 1B: Biofilm production
#################################
biofilmNonhealing = FullData %>% filter(Healed=="Unhealed")
biofilmHealing = FullData %>% filter(Healed=="Healed")
biofilmNonhealingSort = biofilmNonhealing %>% arrange(biofilm) 
biofilmNonhealingSort$IsolateNum = 1:nrow(biofilmNonhealingSort)

biofilmHealingSort = biofilmHealing %>% arrange(biofilm)
biofilmHealingSort$IsolateNum = 1:nrow(biofilmHealingSort)


OrderedBiofilm = rbind(biofilmHealingSort, biofilmNonhealingSort)

biofilm_plot = ggplot(data=OrderedBiofilm, aes(x=Healed, y=IsolateNum, fill=biofilm)) + geom_tile()+scale_fill_gradient(low = "white", high = "orange2")+
  labs(y="Isolates", fill="Biofilm", x="")+theme_classic() +
  theme(axis.title.x = element_text(face="bold",size=15), axis.title.y = element_text(face="bold",size=15), 
        axis.text.y = element_text(face="bold",size=15),axis.text.x = element_text(face="bold",size=15),
        plot.title=element_text(hjust=.5, face="bold", size=20) ) + ggtitle("Biofilm Production") 


# Figure 1C: Siderophore production
#################################
siderophoreNonhealing = FullData %>% filter(Healed=="Unhealed")
siderophoreHealing = FullData %>% filter(Healed=="Healed")

siderophoreNonhealingSort = siderophoreNonhealing %>% arrange(siderophore) 
siderophoreNonhealingSort$IsolateNum = 1:nrow(siderophoreNonhealingSort)
siderophoreHealingSort = siderophoreHealing %>% arrange(siderophore)
siderophoreHealingSort$IsolateNum = 1:nrow(siderophoreHealingSort)


OrderedSiderophore = rbind(siderophoreHealingSort, siderophoreNonhealingSort)

siderophore_plot = ggplot(data=OrderedSiderophore, aes(x=Healed, y=IsolateNum, fill=siderophore)) + geom_tile()+scale_fill_gradient(low = "white", high = "orange2")+
  labs(y="Isolates", fill="Siderophore", x="")+theme_classic() +
  theme(axis.title.x = element_text(face="bold",size=15), axis.title.y = element_text(face="bold",size=15), 
        axis.text.y = element_text(face="bold",size=15),axis.text.x = element_text(face="bold",size=15),
        plot.title=element_text(hjust=.5, face="bold", size=20) ) + ggtitle("Siderophore Production") 



# Figure 1D staphyloxanthin production
######################################

staphyloxanthinNonhealing = FullData %>% filter(Healed=="Unhealed")
staphyloxanthinHealing = FullData %>% filter(Healed=="Healed")
staphyloxanthinNonhealingSort = staphyloxanthinNonhealing %>% arrange(staphyloxanthin) 
staphyloxanthinNonhealingSort$IsolateNum = 1:nrow(staphyloxanthinNonhealingSort)

staphyloxanthinHealingSort = staphyloxanthinHealing %>% arrange(staphyloxanthin)
staphyloxanthinHealingSort$IsolateNum = 1:nrow(staphyloxanthinHealingSort)


Orderedstaphyloxanthin = rbind(staphyloxanthinHealingSort, staphyloxanthinNonhealingSort)

staphyloxanthin_plot = ggplot(data=Orderedstaphyloxanthin, aes(x=Healed, y=IsolateNum, fill=staphyloxanthin)) + geom_tile()+scale_fill_gradient(low = "white", high = "orange2")+
  labs(y="Isolates", fill="Staphyloxanthin", x="")+theme_classic() +
  theme(axis.title.x = element_text(face="bold",size=15), axis.title.y = element_text(face="bold",size=15), 
        axis.text.y = element_text(face="bold",size=15),axis.text.x = element_text(face="bold",size=15),
        plot.title=element_text(hjust=.5, face="bold", size=20) ) + ggtitle("Staphyloxanthin Production") 


# Figure 1E staphylokinase production
######################################

staphylokinaseNonhealing = FullData %>% filter(Healed=="Unhealed")
staphylokinaseHealing = FullData %>% filter(Healed=="Healed")
staphylokinaseNonhealingSort = staphylokinaseNonhealing %>% arrange(staphylokinase) 
staphylokinaseNonhealingSort$IsolateNum = 1:nrow(staphylokinaseNonhealingSort)

staphylokinaseHealingSort = staphylokinaseHealing %>% arrange(staphylokinase)
staphylokinaseHealingSort$IsolateNum = 1:nrow(staphylokinaseHealingSort)

Orderedstaphylokinase = rbind(staphylokinaseHealingSort, staphylokinaseNonhealingSort)

staphylokinase_plot = ggplot(data=Orderedstaphylokinase, aes(x=Healed, y=IsolateNum, fill=staphylokinase)) + geom_tile()+scale_fill_gradient(low = "white", high = "orange2")+
  labs(y="Isolates", fill="Staphylokinase", x="")+theme_classic() +
  theme(axis.title.x = element_text(face="bold",size=15), axis.title.y = element_text(face="bold",size=15), 
        axis.text.y = element_text(face="bold",size=15),axis.text.x = element_text(face="bold",size=15),
        plot.title=element_text(hjust=.5, face="bold", size=20) ) + ggtitle("Staphylokinase Production") 


staphylokinase_plot

pdf(width=24, height= 6.171429, file="Figures/Figure1/Fig1B_E.pdf")
grid.arrange(biofilm_plot, siderophore_plot, staphyloxanthin_plot, staphylokinase_plot, ncol=4)
dev.off()


# Maximum measures per patient of : # logStaphyloxanthin siderophore_normalized logBiofilm staphylokinase_normalized
##############################
Maxed = FullData %>% group_by(patient, Healed, Healed_by_12) %>% summarize(max_staphyloxanthin = max(logStaphyloxanthin), max_biofilm = max(logBiofilm), max_siderophore = max(siderophore_normalized), max_staphylokinase = max(staphylokinase_normalized,na.rm=T))


T_biofilm = t.test(Maxed$max_biofilm~ Maxed$Healed)
p_biofilm = round(T_biofilm$p.value, 4)

T_siderophore = t.test(Maxed$max_siderophore~Maxed$Healed)
p_siderophore = round(T_siderophore$p.value, 4)

T_staphyloxanthin = t.test(Maxed$max_staphyloxanthin~Maxed$Healed)
p_staphyloxanthin = round(T_staphyloxanthin$p.value, 4)

T_staphylokinase = t.test(Maxed$max_staphylokinase~Maxed$Healed)
p_staphylokinase = round(T_staphylokinase$p.value, 4)


# Figure 1F maximum biofilm production per patient
##################################################
biofilm_boxplot = ggplot(data=Maxed, aes(x=Healed, y=max_biofilm)) + geom_boxplot() + labs(x="Healing Outcome", y="Maximum normalized biofilm") + 
  geom_jitter(aes(color=Healed),size=3)+  theme_classic()+
  theme(axis.title.x = element_text(face="bold",size=15),
        axis.title.y = element_text(face="bold",size=15),
        axis.text.y = element_text(face="bold",size=10), 
        axis.text.x = element_text(face="bold",size=10),
        plot.title = element_text(hjust = 0.5, face="bold", size=20))+ scale_color_manual(values=c("grey","orange2")) + annotate(geom="text", label=paste0("p=", p_biofilm), x=2, y=5, size=8)
biofilm_boxplot

# Figure 1G maximum siderophore production per patient
######################################################
siderophore_boxplot = ggplot(data=Maxed, aes(x=Healed, y=max_siderophore)) + geom_boxplot() + labs(x="Healing Outcome", y="Maximum normalized siderophore") + 
  geom_jitter(aes(color=Healed),size=3)+  theme_classic()+
  theme(axis.title.x = element_text(face="bold",size=15),
        axis.title.y = element_text(face="bold",size=15),
        axis.text.y = element_text(face="bold",size=10), 
        axis.text.x = element_text(face="bold",size=10),
        plot.title = element_text(hjust = 0.5, face="bold", size=20))+ scale_color_manual(values=c("grey","orange2")) + annotate(geom="text", label=paste0("p=", p_siderophore), x=2, y=3, size=8)
siderophore_boxplot


# Figure 1H maximum staphyloxanthin production per patient
##########################################################
staphyloxanthin_boxplot = ggplot(data=Maxed, aes(x=Healed, y=max_staphyloxanthin)) + geom_boxplot() + labs(x="Healing Outcome", y="Maximum normalized staphyloxanthin") + 
  geom_jitter(aes(color=Healed),size=3)+  theme_classic()+
  theme(axis.title.x = element_text(face="bold",size=15),
        axis.title.y = element_text(face="bold",size=15),
        axis.text.y = element_text(face="bold",size=10), 
        axis.text.x = element_text(face="bold",size=10),
        plot.title = element_text(hjust = 0.5, face="bold", size=20))+ scale_color_manual(values=c("grey","orange2")) + annotate(geom="text", label=paste0("p=", p_staphyloxanthin), x=2, y=4.75, size=8)
staphyloxanthin_boxplot

# Figure 1I maximum staphylokinase production per patient
##########################################################
staphylokinase_boxplot = ggplot(data=Maxed, aes(x=Healed, y=max_staphylokinase)) + geom_boxplot() + labs(x="Healing Outcome", y="Maximum normalized staphylokinase") + 
  geom_jitter(aes(color=Healed),size=3)+  theme_classic()+
  theme(axis.title.x = element_text(face="bold",size=15),
        axis.title.y = element_text(face="bold",size=15),
        axis.text.y = element_text(face="bold",size=10), 
        axis.text.x = element_text(face="bold",size=10),
        plot.title = element_text(hjust = 0.5, face="bold", size=20))+ scale_color_manual(values=c("grey","orange2")) + annotate(geom="text", label=paste0("p=", p_staphylokinase), x=2, y=3.5, size=8)
staphylokinase_boxplot


pdf(width=24, height= 6.171429, file="Figures/Figure1/Fig1F-I.pdf")
grid.arrange(biofilm_boxplot, siderophore_boxplot, staphyloxanthin_boxplot, staphylokinase_boxplot, ncol=4)
dev.off()

# Figure 1 Legends: LMMs for phenotypes ~ 1|patient + healed 
##############################################################
ZeroOneNormed = FullData
ZeroOneNormed$staphyloxanthinZO = sapply(ZeroOneNormed$staphyloxanthin, function(x) x- min(ZeroOneNormed$staphyloxanthin)) /(max(ZeroOneNormed$staphyloxanthin)-min(ZeroOneNormed$staphyloxanthin))
ZeroOneNormed$staphylokinaseZO = sapply(ZeroOneNormed$staphylokinase, function(x) x- min(ZeroOneNormed$staphylokinase, na.rm=T)) /(max(ZeroOneNormed$staphylokinase, na.rm=T)-min(ZeroOneNormed$staphylokinase, na.rm=T))
ZeroOneNormed$biofilmZO = sapply(ZeroOneNormed$biofilm, function(x) x- min(ZeroOneNormed$biofilm, na.rm=T)) /(max(ZeroOneNormed$biofilm, na.rm=T)-min(ZeroOneNormed$biofilm, na.rm=T))
ZeroOneNormed$siderophoreZO = sapply(ZeroOneNormed$siderophore, function(x) x- min(ZeroOneNormed$siderophore, na.rm=T)) /(max(ZeroOneNormed$siderophore, na.rm=T)-min(ZeroOneNormed$siderophore, na.rm=T))

Xanthin = lme4::lmer( staphyloxanthinZO ~ (1 | patient) + Healed , data=ZeroOneNormed, REML = F)
summary(Xanthin)

siderophoreLMM = lme4::lmer( siderophoreZO ~ (1 | patient) + Healed  , data=ZeroOneNormed, REML = F)
summary(siderophoreLMM)

biofilmLMM = lme4::lmer( biofilmZO ~ (1 | patient) + Healed, data=ZeroOneNormed, REML = F)
summary(biofilmLMM)

staphylokinaseLMM = lme4::lmer( staphylokinaseZO ~ (1 | patient) + Healed , data=ZeroOneNormed, REML = F)
summary(staphylokinaseLMM)


ZeroOneNormed

