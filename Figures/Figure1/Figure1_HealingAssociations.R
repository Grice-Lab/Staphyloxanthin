# Figure 1
# Healing associations
# Modified by Amy Campbell from Amelia McCready-Vangi's and Amy Campbell's scripts
# To reproduce figures 1B-E

library("ggplot2")
library("gridExtra")
library("lme4")
library("dplyr")
library("bestNormalize")

set.seed(19104)
setwd("/Users/amycampbell/Desktop/GriceLabGit/Staphyloxanthin/")

FullData = read.csv("Data/InVitroData/staphyloxanthin_paper_updated.csv")

FullData$Healed = if_else(is.na(FullData$week_healed ), "Unhealed", "Healed")


# For the supplements
FullData$Healed_by_12 = if_else(is.na(FullData$week_healed ) | FullData$week_healed > 12, "Unhealed at 12 Weeks", "Healed at 12 Weeks")

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

# Figure 1B(i): Biofilm production
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

biofilm_plot_viridis = ggplot(data=OrderedBiofilm, aes(x=Healed, y=IsolateNum, fill=biofilm)) + geom_tile()+scale_fill_viridis(option="plasma")+
  labs(y="Isolates", fill="Biofilm", x="")+theme_classic() +
  theme(axis.title.x = element_text(face="bold",size=15), axis.title.y = element_text(face="bold",size=15), 
        axis.text.y = element_text(face="bold",size=15),axis.text.x = element_text(face="bold",size=15),
        plot.title=element_text(hjust=.5, face="bold", size=20) ) + ggtitle("Biofilm Production") 

# Figure 1C(i): Siderophore production
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



# Figure 1D(i) staphylokinase production
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



pdf(width=24, height= 6.171429, file="Figures/Figure1/Fig1B_E.pdf")
grid.arrange(biofilm_plot, siderophore_plot, staphyloxanthin_plot, staphylokinase_plot, ncol=4)
dev.off()




# Figure 1E(i) staphyloxanthin production
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


# Maximum measures per patient of : # logStaphyloxanthin siderophore_normalized logBiofilm staphylokinase_normalized
######################################################################################################################
Maxed = FullData %>% group_by(patient, Healed, Healed_by_12) %>% summarize(max_staphyloxanthin = max(logStaphyloxanthin), max_biofilm = max(logBiofilm), max_siderophore = max(siderophore_normalized), max_staphylokinase = max(staphylokinase_normalized,na.rm=T))


T_biofilm = t.test(Maxed$max_biofilm~ Maxed$Healed)
p_biofilm = round(T_biofilm$p.value, 4)

T_siderophore = t.test(Maxed$max_siderophore~Maxed$Healed)
p_siderophore = round(T_siderophore$p.value, 4)

T_staphyloxanthin = t.test(Maxed$max_staphyloxanthin~Maxed$Healed)
p_staphyloxanthin = round(T_staphyloxanthin$p.value, 4)

T_staphylokinase = t.test(Maxed$max_staphylokinase~Maxed$Healed)
p_staphylokinase = round(T_staphylokinase$p.value, 4)


# Figure 1B(ii) maximum biofilm production per patient
##################################################
biofilm_boxplot = ggplot(data=Maxed, aes(x=Healed, y=max_biofilm)) + geom_boxplot() + labs(x="Healing Outcome", y="Maximum normalized biofilm") + 
  geom_jitter(aes(color=Healed),size=3, height=0)+  theme_classic()+
  theme(axis.title.x = element_text(face="bold",size=15),
        axis.title.y = element_text(face="bold",size=15),
        axis.text.y = element_text(face="bold",size=10), 
        axis.text.x = element_text(face="bold",size=10),
        plot.title = element_text(hjust = 0.5, face="bold", size=20))+ scale_color_manual(values=c("grey","orange2")) + annotate(geom="text", label=paste0("p=", p_biofilm), x=2, y=5, size=8)
biofilm_boxplot

# Figure 1C(ii) maximum siderophore production per patient
######################################################
siderophore_boxplot = ggplot(data=Maxed, aes(x=Healed, y=max_siderophore)) + geom_boxplot() + labs(x="Healing Outcome", y="Maximum normalized siderophore") + 
  geom_jitter(aes(color=Healed),size=3,height=0)+  theme_classic()+
  theme(axis.title.x = element_text(face="bold",size=15),
        axis.title.y = element_text(face="bold",size=15),
        axis.text.y = element_text(face="bold",size=10), 
        axis.text.x = element_text(face="bold",size=10),
        plot.title = element_text(hjust = 0.5, face="bold", size=20))+ scale_color_manual(values=c("grey","orange2")) + annotate(geom="text", label=paste0("p=", p_siderophore), x=2, y=3, size=8)
siderophore_boxplot

# Figure 1d(ii) maximum staphylokinase production per patient
##########################################################
staphylokinase_boxplot = ggplot(data=Maxed, aes(x=Healed, y=max_staphylokinase)) + geom_boxplot() + labs(x="Healing Outcome", y="Maximum normalized staphylokinase") + 
  geom_jitter(aes(color=Healed),size=3, height=0)+  theme_classic()+
  theme(axis.title.x = element_text(face="bold",size=15),
        axis.title.y = element_text(face="bold",size=15),
        axis.text.y = element_text(face="bold",size=10), 
        axis.text.x = element_text(face="bold",size=10),
        plot.title = element_text(hjust = 0.5, face="bold", size=20))+ scale_color_manual(values=c("grey","orange2")) + annotate(geom="text", label=paste0("p=", p_staphylokinase), x=2, y=3.5, size=8)
staphylokinase_boxplot

# Figure 1e(ii) maximum staphyloxanthin production per patient
##########################################################
staphyloxanthin_boxplot = ggplot(data=Maxed, aes(x=Healed, y=max_staphyloxanthin)) + geom_boxplot() + labs(x="Healing Outcome", y="Maximum normalized staphyloxanthin") + 
  geom_jitter(aes(color=Healed),size=3,height=0)+  theme_classic()+
  theme(axis.title.x = element_text(face="bold",size=15),
        axis.title.y = element_text(face="bold",size=15),
        axis.text.y = element_text(face="bold",size=10), 
        axis.text.x = element_text(face="bold",size=10),
        plot.title = element_text(hjust = 0.5, face="bold", size=20))+ scale_color_manual(values=c("grey","orange2")) + annotate(geom="text", label=paste0("p=", p_staphyloxanthin), x=2, y=4.75, size=8)
staphyloxanthin_boxplot



pdf(width=24, height= 6.171429, file="Figures/Figure1/Fig1F-I.pdf")
grid.arrange(biofilm_boxplot, siderophore_boxplot, staphyloxanthin_boxplot, staphylokinase_boxplot, ncol=4)
dev.off()

# Figure S1a max biofilm production per patient -- 12 weeks status
##################################################################
T_biofilm12 = t.test(Maxed$max_biofilm~ Maxed$Healed_by_12)
p_biofilm12 = round(T_biofilm12$p.value, 4)

T_siderophore12 = t.test(Maxed$max_siderophore~Maxed$Healed_by_12)
p_siderophore12 = round(T_siderophore12$p.value, 4)


T_staphylokinase12 = t.test(Maxed$max_staphylokinase~Maxed$Healed_by_12)
p_staphylokinase12 = round(T_staphylokinase12$p.value, 4)

T_staphyloxanthin12 = t.test(Maxed$max_staphyloxanthin~Maxed$Healed_by_12)
p_staphyloxanthin12 = round(T_staphyloxanthin12$p.value, 4)

# Figure S1A maximum biofilm production per patient
##################################################
biofilm_boxplotS1 = ggplot(data=Maxed, aes(x=Healed_by_12, y=max_biofilm)) + geom_boxplot() + labs( y="Maximum normalized biofilm") + 
  geom_jitter(aes(color=Healed_by_12),size=3,height=0)+  theme_classic()+
  theme(axis.title.y = element_text(face="bold",size=15),
        axis.text.y = element_text(face="bold",size=10), 
        axis.text.x = element_text(face="bold",size=10),
        axis.title.x=element_blank(),
        plot.title = element_text(hjust = 0.5, face="bold", size=20))+ scale_color_manual(values=c("grey","orange2")) + annotate(geom="text", label=paste0("p=", p_biofilm12), x=2, y=5, size=8)
biofilm_boxplotS1

# Figure 1C(ii) maximum siderophore production per patient
######################################################
siderophore_boxplotS1 = ggplot(data=Maxed, aes(x=Healed_by_12, y=max_siderophore)) + geom_boxplot() + labs( y="Maximum normalized siderophore") + 
  geom_jitter(aes(color=Healed_by_12),size=3, height=0)+  theme_classic()+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(face="bold",size=15),
        axis.text.y = element_text(face="bold",size=10), 
        axis.text.x = element_text(face="bold",size=10),
        plot.title = element_text(hjust = 0.5, face="bold", size=20))+ scale_color_manual(values=c("grey","orange2")) + annotate(geom="text", label=paste0("p=", p_siderophore12), x=2, y=3, size=8)
siderophore_boxplotS1

# Figure 1d(ii) maximum staphylokinase production per patient
##########################################################
staphylokinase_boxplotS1 = ggplot(data=Maxed, aes(x=Healed_by_12, y=max_staphylokinase)) + geom_boxplot() + labs( y="Maximum normalized staphylokinase") + 
  geom_jitter(aes(color=Healed_by_12),size=3, height=0)+  theme_classic()+
  theme(axis.title.y = element_text(face="bold",size=15),
        axis.text.y = element_text(face="bold",size=10), 
        axis.text.x = element_text(face="bold",size=10),
        axis.title.x=element_blank(),
        plot.title = element_text(hjust = 0.5, face="bold", size=20))+ scale_color_manual(values=c("grey","orange2")) + annotate(geom="text", label=paste0("p=", p_staphylokinase12), x=2, y=3.5, size=8)
staphylokinase_boxplotS1

# Figure 1e(ii) maximum staphyloxanthin production per patient
##########################################################
staphyloxanthin_boxplotS1 = ggplot(data=Maxed, aes(x=Healed_by_12, y=max_staphyloxanthin)) + geom_boxplot() + labs( y="Maximum normalized staphyloxanthin") + 
  geom_jitter(aes(color=Healed_by_12),size=3, height=0)+  theme_classic()+
  theme(axis.title.y = element_text(face="bold",size=15),
        axis.text.y = element_text(face="bold",size=10), 
        axis.title.x=element_blank(),
        axis.text.x = element_text(face="bold",size=10),
        plot.title = element_text(hjust = 0.5, face="bold", size=20))+ scale_color_manual(values=c("grey","orange2")) + annotate(geom="text", label=paste0("p=", p_staphyloxanthin12), x=2, y=4.75, size=8)
staphyloxanthin_boxplotS1

pdf(width=11, height= 12, file="Figures/Figure1/FigS1A-D.pdf")
grid.arrange(biofilm_boxplotS1, siderophore_boxplotS1, staphylokinase_boxplotS1, staphyloxanthin_boxplotS1, ncol=2)
dev.off()



# Table S2: LMMs for phenotypes ~ 1|patient + healed 
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



# Figure S_? : Repeat de-duplication but use MEDIAN values instead to mitigate effect of extreme values 
#######################################################################################################
Medians = FullData %>% group_by(patient, Healed, Healed_by_12) %>% summarize(med_staphyloxanthin = median(logStaphyloxanthin), med_biofilm = median(logBiofilm), med_siderophore = median(siderophore_normalized), med_staphylokinase = median(staphylokinase_normalized,na.rm=T))

T_median_biofilm = t.test(Medians$med_biofilm~ Medians$Healed)
p_median_biofilm = round(T_median_biofilm$p.value, 4)

T_median_siderophore = t.test(Medians$med_siderophore~Medians$Healed)
p_median_siderophore = round(T_median_siderophore$p.value, 4)

T_median_staphyloxanthin = t.test(Medians$med_staphyloxanthin~Medians$Healed)
p_median_staphyloxanthin = round(T_median_staphyloxanthin$p.value, 4)

T_median_staphylokinase = t.test(Medians$med_staphylokinase~Medians$Healed)
p_median_staphylokinase = round(T_median_staphylokinase$p.value, 4)



# Figure d(ii) median biofilm production per patient
##################################################
biofilm_boxplot_median = ggplot(data=Medians, aes(x=Healed, y=med_biofilm)) + geom_boxplot() + labs(x="Healing Outcome", y="Median normalized biofilm") + 
  geom_jitter(aes(color=Healed),size=3, height=0)+  theme_classic()+
  theme(axis.title.x = element_text(face="bold",size=15),
        axis.title.y = element_text(face="bold",size=15),
        axis.text.y = element_text(face="bold",size=10), 
        axis.text.x = element_text(face="bold",size=10),
        plot.title = element_text(hjust = 0.5, face="bold", size=20))+ scale_color_manual(values=c("grey","orange2")) + annotate(geom="text", label=paste0("p=", p_median_biofilm), x=2, y=5, size=8)
biofilm_boxplot_median

# Figure d(ii) median siderophore production per patient
######################################################
siderophore_boxplot_median = ggplot(data=Medians, aes(x=Healed, y=med_siderophore)) + geom_boxplot() + labs(x="Healing Outcome", y="Median normalized siderophore") + 
  geom_jitter(aes(color=Healed),size=3,height=0)+  theme_classic()+
  theme(axis.title.x = element_text(face="bold",size=15),
        axis.title.y = element_text(face="bold",size=15),
        axis.text.y = element_text(face="bold",size=10), 
        axis.text.x = element_text(face="bold",size=10),
        plot.title = element_text(hjust = 0.5, face="bold", size=20))+ scale_color_manual(values=c("grey","orange2")) + annotate(geom="text", label=paste0("p=", p_median_siderophore), x=2, y=3, size=8)
siderophore_boxplot_median

# Figure d(ii) median staphylokinase production per patient
##########################################################
staphylokinase_boxplot_median = ggplot(data=Medians, aes(x=Healed, y=med_staphylokinase)) + geom_boxplot() + labs(x="Healing Outcome", y="Median normalized staphylokinase") + 
  geom_jitter(aes(color=Healed),size=3, height=0)+  theme_classic()+
  theme(axis.title.x = element_text(face="bold",size=15),
        axis.title.y = element_text(face="bold",size=15),
        axis.text.y = element_text(face="bold",size=10), 
        axis.text.x = element_text(face="bold",size=10),
        plot.title = element_text(hjust = 0.5, face="bold", size=20))+ scale_color_manual(values=c("grey","orange2")) + annotate(geom="text", label=paste0("p=", p_median_staphylokinase), x=2, y=3.5, size=8)
staphylokinase_boxplot

# Figure e(ii) median staphyloxanthin production per patient
##########################################################
staphyloxanthin_boxplot_median = ggplot(data=Medians, aes(x=Healed, y=med_staphyloxanthin)) + geom_boxplot() + labs(x="Healing Outcome", y="Median normalized staphyloxanthin") + 
  geom_jitter(aes(color=Healed),size=3, height=0)+  theme_classic()+
  theme(axis.title.x = element_text(face="bold",size=15),
        axis.title.y = element_text(face="bold",size=15),
        axis.text.y = element_text(face="bold",size=10), 
        axis.text.x = element_text(face="bold",size=10),
        plot.title = element_text(hjust = 0.5, face="bold", size=20))+ scale_color_manual(values=c("grey","orange2")) + annotate(geom="text", label=paste0("p=", p_median_staphyloxanthin), x=2, y=4.75, size=8)

pdf(width=9, height= 12, file="Figures/Figure1/FigS1A-D_Medians.pdf")
grid.arrange(biofilm_boxplot_median, siderophore_boxplot_median, staphylokinase_boxplot_median, staphyloxanthin_boxplot_median, ncol=2)
dev.off()



# Figure S_? : Repeat de-duplication but use MEAN values instead to mitigate effect of extreme values 
#######################################################################################################
Means = FullData %>% group_by(patient, Healed, Healed_by_12) %>% summarize(mean_staphyloxanthin = mean(logStaphyloxanthin, na.rm=T), mean_biofilm = mean(logBiofilm, na.rm=T), mean_siderophore = mean(siderophore_normalized,na.rm=T), mean_staphylokinase = mean(staphylokinase_normalized,na.rm=T))

T_mean_biofilm = t.test(Means$mean_biofilm~ Means$Healed)
p_mean_biofilm = round(T_mean_biofilm$p.value, 4)

T_mean_siderophore = t.test(Means$mean_siderophore~Means$Healed)
p_mean_siderophore = round(T_mean_siderophore$p.value, 4)

T_mean_staphyloxanthin = t.test(Means$mean_staphyloxanthin~Means$Healed)
p_mean_staphyloxanthin = round(T_mean_staphyloxanthin$p.value, 4)

T_mean_staphylokinase = t.test(Means$mean_staphylokinase~Means$Healed)
p_mean_staphylokinase = round(T_mean_staphylokinase$p.value, 4)


# Figure d(ii) mean biofilm production per patient
##################################################

biofilm_boxplot_mean = ggplot(data=Means, aes(x=Healed, y=mean_biofilm)) + geom_boxplot() + labs(x="Healing Outcome", y="Mean normalized biofilm") + 
  geom_jitter(aes(color=Healed),size=3, height=0)+  theme_classic()+
  theme(axis.title.x = element_text(face="bold",size=15),
        axis.title.y = element_text(face="bold",size=15),
        axis.text.y = element_text(face="bold",size=10), 
        axis.text.x = element_text(face="bold",size=10),
        plot.title = element_text(hjust = 0.5, face="bold", size=20))+ scale_color_manual(values=c("grey","orange2")) + annotate(geom="text", label=paste0("p=", p_mean_biofilm), x=2, y=5, size=8)
biofilm_boxplot_mean

# Figure d(ii) mean siderophore production per patient
######################################################
siderophore_boxplot_mean = ggplot(data=Means, aes(x=Healed, y=mean_siderophore)) + geom_boxplot() + labs(x="Healing Outcome", y="Mean normalized siderophore") + 
  geom_jitter(aes(color=Healed),size=3, height=0)+  theme_classic()+
  theme(axis.title.x = element_text(face="bold",size=15),
        axis.title.y = element_text(face="bold",size=15),
        axis.text.y = element_text(face="bold",size=10), 
        axis.text.x = element_text(face="bold",size=10),
        plot.title = element_text(hjust = 0.5, face="bold", size=20))+ scale_color_manual(values=c("grey","orange2")) + annotate(geom="text", label=paste0("p=", p_mean_siderophore), x=2, y=3, size=8)
siderophore_boxplot_mean

# Figure d(ii) mean staphylokinase production per patient
##########################################################
staphylokinase_boxplot_mean = ggplot(data=Means, aes(x=Healed, y=mean_staphylokinase)) + geom_boxplot() + labs(x="Healing Outcome", y="Mean normalized staphylokinase") + 
  geom_jitter(aes(color=Healed),size=3, height=0)+  theme_classic()+
  theme(axis.title.x = element_text(face="bold",size=15),
        axis.title.y = element_text(face="bold",size=15),
        axis.text.y = element_text(face="bold",size=10), 
        axis.text.x = element_text(face="bold",size=10),
        plot.title = element_text(hjust = 0.5, face="bold", size=20))+ scale_color_manual(values=c("grey","orange2")) + annotate(geom="text", label=paste0("p=", p_mean_staphylokinase), x=2, y=3.5, size=8)
staphylokinase_boxplot

# Figure e(ii) mean staphyloxanthin production per patient
##########################################################
staphyloxanthin_boxplot_mean = ggplot(data=Means, aes(x=Healed, y=mean_staphyloxanthin)) + geom_boxplot() + labs(x="Healing Outcome", y="Mean normalized staphyloxanthin") + 
  geom_jitter(aes(color=Healed),size=3, height=0)+  theme_classic()+
  theme(axis.title.x = element_text(face="bold",size=15),
        axis.title.y = element_text(face="bold",size=15),
        axis.text.y = element_text(face="bold",size=10), 
        axis.text.x = element_text(face="bold",size=10),
        plot.title = element_text(hjust = 0.5, face="bold", size=20))+ scale_color_manual(values=c("grey","orange2")) + annotate(geom="text", label=paste0("p=", p_mean_staphyloxanthin), x=2, y=4.75, size=8) 
staphyloxanthin_boxplot_mean


pdf(width=9, height= 12, file="Figures/Figure1/FigS1A-D_Means.pdf")
grid.arrange(biofilm_boxplot_mean, siderophore_boxplot_mean, staphylokinase_boxplot, staphyloxanthin_boxplot_mean, ncol=2)
dev.off()


# What happens when we consider the first isolates collected 
# post-debridement? (in the time post-debridement, wounds are
# supposed to be in the inflammatory/acute healing phase)
# (only includes 45 patients because 15 of them only had
# isolates from timepoint 0, which is pre-debridement)
#########################################################

minVisitDF = data.frame()
maxVisitDF = data.frame()

# Aside from the pre-debridement timepoint, 
# Take first timepoint there's a S. aureus isolate from in each patient
# if there are multiple timepoints, average their values for the phenotypes together
NonZero = FullData %>% filter(visit!=0)
for(p in unique(FullData$patient)){
  PatientData = NonZero %>% subset(patient==p)
  minvis = min(PatientData$visit)
  
  maxvis = max(PatientData$visit)

  
  PatientDataMinimumVis = PatientData %>% filter(visit==minvis)
  PatientDataMaximumVis = PatientData %>% filter(visit==maxvis)
  
  stx=mean(PatientDataMinimumVis$logStaphyloxanthin, na.rm=T)
  stx_maxvisit=mean(PatientDataMaximumVis$logStaphyloxanthin, na.rm=T)
  
  staphylokinase=mean(PatientDataMinimumVis$staphylokinase_normalized, na.rm=T)
  staphylokinase_maxvisit = mean(PatientDataMaximumVis$logStaphyloxanthin, na.rm=T)
  
  biofilm=mean(PatientDataMinimumVis$logBiofilm, na.rm=T)
  biofilm_maxvisit =  mean(PatientDataMaximumVis$logBiofilm, na.rm=T)
  
  siderophores=mean(PatientDataMinimumVis$siderophore_normalized, na.rm=T)
  siderophores_maxvisit = mean(PatientDataMaximumVis$siderophore_normalized, na.rm=T)
  
  rowitem = c(p, minvis, stx, staphylokinase, biofilm, siderophores )
  rowitem_max = c(p, maxvis, stx_maxvisit, staphylokinase_maxvisit, biofilm_maxvisit, siderophores_maxvisit )
  minVisitDF = rbind(minVisitDF, rowitem )
  maxVisitDF = rbind(maxVisitDF, rowitem_max )
}

colnames(minVisitDF) = c("patient", "Visit", "Staphyloxanthin", "Staphylokinase", "Biofilm", "Siderophores")

minVisitDF = minVisitDF %>% left_join(FullData %>% select(patient,Healed ), by="patient") %>% unique()

minVisitDF = minVisitDF %>% filter(!is.na(Staphyloxanthin))

T_firstIsolate_staphyloxanthin = t.test(minVisitDF$Staphyloxanthin~minVisitDF$Healed)
p_firstIsolate_staphyloxanthin = round(T_firstIsolate_staphyloxanthin$p.value, 4)

T_firstIsolate_Staphylokinase = t.test(minVisitDF$Staphylokinase~minVisitDF$Healed)
p_firstIsolate_Staphylokinase = round(T_firstIsolate_Staphylokinase$p.value, 4)

T_firstIsolate_biofilm = t.test(minVisitDF$Biofilm~minVisitDF$Healed)
p_firstIsolate_biofilm = round(T_firstIsolate_biofilm$p.value, 4)

T_firstIsolate_siderophore = t.test(minVisitDF$Siderophores~minVisitDF$Healed)
p_firstIsolate_siderophore = round(T_firstIsolate_siderophore$p.value, 4)



MeanByVisit = FullData %>% group_by(patient,visit) %>% summarize(patient=patient, visit=visit,
                                                                 staphyloxanthinmean=mean(logStaphyloxanthin),
                                                                 biofilmmean=mean(logBiofilm),
                                                                 siderophoremean=mean(siderophore_normalized), 
                                                                 staphylokinasemean=mean(staphylokinase_normalized, na.rm=T),
                                                                 week_healed=week_healed)

MeanByVisit = MeanByVisit %>% unique()
MeanByVisit$week = MeanByVisit$week_healed

WeekPalette =RColorBrewer::brewer.pal(11, "Spectral")

darkpurple =RColorBrewer::brewer.pal(11, "PuOr")[10]
WeekPalette = c("#8B0000", WeekPalette)

WeekPalette = c("black", WeekPalette)
WeekPalette = append(WeekPalette, darkpurple)
colormap=data.frame(colors=WeekPalette, week=c(0,2,4,6,8,10,12,14,16,18,20,22,24,26))


MeanByVisit_healed = MeanByVisit %>% filter(!is.na(week_healed))
MeanByVisit_unhealed = MeanByVisit %>% filter(is.na(week_healed))


summary(lm(MeanByVisit_unhealed$staphyloxanthinmean~ MeanByVisit_unhealed$visit))
summary(lm(MeanByVisit_unhealed$biofilmmean~ MeanByVisit_unhealed$visit))
summary(lm(MeanByVisit_unhealed$staphylokinasemean~ MeanByVisit_unhealed$visit))
summary(lm(MeanByVisit_unhealed$siderophoremean~ MeanByVisit_unhealed$visit))

summary(lm(MeanByVisit_healed$staphyloxanthinmean~ MeanByVisit_healed$visit))
summary(lm(MeanByVisit_healed$biofilmmean~ MeanByVisit_healed$visit))
summary(lm(MeanByVisit_healed$staphylokinasemean~ MeanByVisit_healed$visit))
summary(lm(MeanByVisit_healed$siderophoremean~ MeanByVisit_healed$visit))

summary(lm(MeanByVisit$staphyloxanthinmean~ MeanByVisit$visit))


MeanByVisit$Healed = if_else(is.na(MeanByVisit$week_healed), "Unhealed", "Healed")
ggplot(MeanByVisit, aes(x=visit, y=staphyloxanthinmean, color=factor(Healed))) + geom_point() + scale_color_manual(values=c("dodgerblue", "darkorange")) + theme_classic()


