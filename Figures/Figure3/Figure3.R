# Figure 3
# adapted from ARM's scripts

library("survminer")
library("ggplot2")
library("dplyr")

Palette5 = (RColorBrewer::brewer.pal(10, "Spectral"))[c(4, 1, 9, 8, 10)]
Palette3 = (RColorBrewer::brewer.pal(10, "Spectral"))[c(9, 8,10)]


woundsize= read.csv("Data/InVivoData/woundsizes.csv")
PilotExp = read.csv("Data/InVivoData/pilot_wound_measurements.csv")
kaplan<-read.csv("Data/InVivoData/mouse_kaplan_05.02.22.csv")
CFU<-read.csv("Data/InVivoData/CFU_data.csv")

wound_order <- c("PBS","LAC", "LAC_dCrtN_1", "SA925", "SA1088")

WilcoxResults = pairwise.wilcox.test(woundsize$day_14_change, woundsize$group, p.adjust.method = "none", exact=F)
WilcoxResultsPValues = WilcoxResults$p.value

pvalPBS_LAC = "6.165e-07"
pvalPBS_CrtN = "4.979e-05"
pvalPBS_SA925 = "1.263e-05"
pvalPBS_SA1088 = round(WilcoxResultsPValues["SA1088", "PBS"], 4)

pvalLAC_DcrtN = round(WilcoxResultsPValues["LAC_dCrtN_1", "LAC"],4)
pvalSA925_SA1088 = "7.4e-04"

Woundplot = ggplot(data=woundsize, aes(x=group, y=day_14_change))+geom_boxplot(size=.25, outlier.shape = NA)+
  geom_jitter(aes(color=group, shape=sex),size=2, width=.2)+
  labs(x="group", y="% original wound size at day 14")+ theme_classic()+
  theme(axis.title.x = element_text(face="bold",size=15))+
  theme(axis.title.y = element_text(face="bold",size=15))+
  theme(axis.text.y = element_text(face="bold",size=15))+
  theme(axis.text.x = element_text(face="bold",size=15))+
  theme(plot.title = element_text(hjust = 0.5, face="bold", size=20))+
  scale_color_manual(values=(RColorBrewer::brewer.pal(10, "Spectral"))[c(4, 1, 9, 8, 10)])+
  scale_x_discrete(limits = wound_order) +
  annotate(geom="text", x=1.5, y=365, label=paste0("p=", pvalPBS_LAC), size=5) + 
  annotate(geom="text", x=2, y=348, label=paste0("p=", pvalPBS_CrtN), size=5) + 
  annotate(geom="text", x=2.5, y=331, label=paste0("p=", pvalPBS_SA925), size=5) + 
  annotate(geom="text", x=3, y=314, label=paste0("p=", pvalPBS_SA1088), size=5) + 
  annotate(geom="text", x=2.5, y=297, label=paste0("p=", pvalLAC_DcrtN), size=5) + 
  annotate(geom="text", x=4.5, y=200, label=paste0("p=", pvalSA925_SA1088),size=5)

ggsave(Woundplot, file="Figures/Figure3/Figure3b.pdf", width=9.6, height=5.5)



# Potential figure 3C: mean % size of original by day
#####################################################
woundsize$ID = paste(woundsize$mouse, woundsize$left.or.right, woundsize$sex, sep="_")

woundsize$day0change = 100

Day0s = woundsize %>% select(ID, mouse, sex, left.or.right, group, day0change)
colnames(Day0s) = c("ID","mouse", "sex", "left.or.right", "group", "PctOriginal")
Day0s$Day = 0

Day3s = woundsize %>% select(ID,mouse, sex, left.or.right, group,day_3_change)
colnames(Day3s) = c("ID", "mouse","sex", "left.or.right", "group", "PctOriginal")
Day3s$Day = 3

BigDF = rbind(Day0s, Day3s)

Day7s = woundsize %>% select(ID, mouse, sex, left.or.right, group,day_7_change)
colnames(Day7s) = c("ID", "mouse", "sex", "left.or.right", "group", "PctOriginal")
Day7s$Day = 7

BigDF = rbind(BigDF, Day7s)

Day14s = woundsize %>% select(ID,mouse, sex, left.or.right, group,day_14_change)
colnames(Day14s) = c("ID","mouse", "sex", "left.or.right", "group", "PctOriginal")
Day14s$Day = 14
BigDF = rbind(BigDF, Day14s)

BigDF = BigDF %>% group_by(group, Day) %>% summarize(Day=mean(Day), MeanPct = mean(PctOriginal,na.rm=T ), sdPct = sd(PctOriginal, na.rm=T) )

AllTimePoints = ggplot(BigDF, aes(x=Day, y=MeanPct, color=group, group=group)) +  geom_point(position=position_dodge(.5)) + geom_line(position=position_dodge(.5)) +
  geom_errorbar(aes(ymin=MeanPct-sdPct, ymax=MeanPct+sdPct,),position=position_dodge(.5)) + scale_color_manual(values=(RColorBrewer::brewer.pal(10, "Spectral"))[c(4, 1, 9, 8, 10)]) + 
  scale_x_continuous(breaks=c(0,3,7,14)) + theme_classic() + 
  theme(axis.title.x = element_text(face="bold",size=15),
        axis.title.y = element_text(face="bold",size=15),
        axis.text.y = element_text(face="bold",size=15),
        axis.text.x = element_text(face="bold",size=15)) + labs(x="Day", y="% Original Wound Size")


ggsave(AllTimePoints, file="Figures/Figure3/Figure3C.pdf", width=9.6, height=5.5)



# Figure S3a: Pilot experiment wound healing results
####################################################

PilotPairwiseWilcox = pairwise.wilcox.test(PilotExp$change_14, PilotExp$group, pool.sd=F, p.adjust.method = "none", exact=F)

SA1088_PBS = round((PilotPairwiseWilcox$p.value)["SA1088", "PBS"],4)
SA925_PBS = round((PilotPairwiseWilcox$p.value)["SA925", "PBS"],4)
SA925_SA1088 = round((PilotPairwiseWilcox$p.value)["SA925", "SA1088"],4)


PilotPlot = ggplot(PilotExp, aes(x=group, y=change_14)) + geom_boxplot(size=.25) + geom_jitter(aes(color=group),size=3)+theme_classic()+
  theme(axis.title.y = element_text(face="bold",size=15),
        axis.title.x = element_text(face="bold",size=15),
        axis.text.y = element_text(face="bold",size=15),
        axis.text.x = element_text(face="bold",size=15))  + scale_color_manual(values=Palette3) +
  annotate(geom="text", x=1.5, y=100, label=paste0("p=", SA1088_PBS), size=5) +
  annotate(geom="text", x=2, y=420, label=paste0("p=", SA925_PBS), size=5) +
  annotate(geom="text", x=2.5, y=445, label=paste0("p=", SA925_SA1088), size=5) + labs(x="Strain", y="% original wound size at day 14")

ggsave(PilotPlot, file="Figures/Figure3/FigureS3a.pdf", width=7, height=5)


# Figure S3b. Survival probabilities 
####################################

model_fit <- survfit(Surv(day, number) ~ group, data = kaplan)

SurvivalCurve = ggsurvplot(model_fit, data = kaplan, censor.shape="|",
           censor.size = 4, conf.int = TRUE, conf.int.style = "step",size=.5,
           xlab = "Time in days")

pdf(file="Figures/Figure3/FigureS3b.pdf", width=7, height=6)
SurvivalCurve
dev.off()


# Figure S3c. Total CFUs 
colnames(CFU)<-c("mouse","group","average","weight","cfu_gram","SA_gram")
CFU

wilcoxCFUTotal = pairwise.wilcox.test(log10(CFU$cfu_gram), CFU$group, pool.sd=F, p.adjust.method = "none", exact=F)
USA300pvaluesCFUTotal = round((wilcoxCFUTotal$p.value)["USA300 LAC", "USA300 dCrtN"], 4)
ClinicalpvaluesCFUTotal = round((wilcoxCFUTotal$p.value)["SA925", "SA1088"], 4)
PBS_USA300 =  round((wilcoxCFUTotal$p.value)["USA300 LAC", "USA300 dCrtN"], 4)
TotalCFUPlot = ggplot(data=CFU, aes(x=group, y=cfu_gram))+geom_boxplot(size=.25)+
  geom_jitter(aes(color=group),size=2)+
  labs(y="Total CFUs per gram")+
  theme_classic()+
  theme(axis.title.x = element_text(face="bold",size=15))+
  theme(axis.title.y = element_text(face="bold",size=15))+
  theme(axis.text.y = element_text(face="bold",size=15))+
  theme(axis.text.x = element_text(face="bold",size=15))+
  theme(plot.title = element_text(hjust = 0.5, face="bold", size=20))+
  scale_color_manual(values=Palette5)+
  scale_y_continuous(trans='log10')+scale_x_discrete(limits = c("PBS", "USA300 LAC", "USA300 dCrtN", "SA925", "SA1088")) +
  annotate(geom="text", x= 2.5, y=1e07, label=paste0("p=",USA300pvaluesCFUTotal),size=5)+
  annotate(geom="text", x= 4.5, y=3e07, label=paste0("p=",ClinicalpvaluesCFUTotal),size=5) 
  

ggsave(TotalCFUPlot, file="Figures/Figure3/FigureS3c.pdf", width=8, height=5)

# Figure S3d: CFUs of S. aureus at endpoint of wounding experiment
###################################################################
wilcoxCFUaureus = pairwise.wilcox.test(CFU$SA_gram, CFU$group, pool.sd=F, p.adjust.method = "none", exact=F)
USA300pvaluesCFUaureus = round((wilcoxCFUaureus$p.value)["USA300 LAC", "USA300 dCrtN"], 4)
ClinicalpvaluesCFUaureus = round((wilcoxCFUaureus$p.value)["SA925", "SA1088"], 4)



AureusCFUPlot = ggplot(data=CFU, aes(x=group, y=SA_gram+.00001))+geom_boxplot(size=.25)+
  geom_jitter(aes(color=group),size=2)+
  labs(y="S. aureus CFUs per gram")+
  theme_classic()+
  theme(axis.title.x = element_text(face="bold",size=15))+
  theme(axis.title.y = element_text(face="bold",size=15))+
  theme(axis.text.y = element_text(face="bold",size=15))+
  theme(axis.text.x = element_text(face="bold",size=15))+
  theme(plot.title = element_text(hjust = 0.5, face="bold", size=20))+
  scale_color_manual(values=Palette5)+
  scale_y_continuous(trans='log10')+scale_x_discrete(limits = c("PBS", "USA300 LAC", "USA300 dCrtN", "SA925", "SA1088")) +
  annotate(geom="text", x= 2.5, y=2e07, label=paste0("p=",USA300pvaluesCFUaureus),size=5)+
  annotate(geom="text", x= 4.5, y=4e07, label=paste0("p=",ClinicalpvaluesCFUaureus),size=5) 
ggsave(AureusCFUPlot, file="Figures/Figure3/FigureS3d.pdf", width=8, height=5)



