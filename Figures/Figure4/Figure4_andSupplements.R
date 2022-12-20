# Amy Campbell
# Figure 4a-b & supplement
## Figures 4a-b, supplements
# Modified from Victoria Lovins' original scripts analyzing these data in late 2021
# Experiments 8(black-6 mice) and 11(db/db mice) performed by VL 

library(tidyverse)
library(reshape2)
library(ggplot2)
library(ggpubr)

FlowData = "Data/InVivoData/Neutrophil_Flow_Experiments_All.csv"
flowdata = read.csv(FlowData)

neutrophils<-read.csv("Data/InVitroData/neutrophil_survival.csv")
blood<-read.csv("Data/InVitroData/whole_blood_survival.csv")


# Figure 4a: whole blood survival of S. aureus 
##############################################


blood_order<-c("USA300", "USA300 dCrtN","SA925","SA1088")
TTestBlood = pairwise.t.test(blood$survival, blood$strain, pool.sd=F, p.adjust.method = "none")
pvalsBlood = TTestBlood$p.value


BloodSurvivalPlot = ggplot(blood, aes(x=strain,y=survival))+geom_boxplot(fill="dodgerblue", width = 0.5, size=.25)+
  geom_jitter(width=.15)+
  labs(y="% Survival After Whole Blood",x="strain")+
  theme(plot.title = element_text(hjust = 0.5, face="bold", size=20),
  axis.title.x = element_text(face="bold",size=15),
  axis.title.y = element_text(face="bold",size=15),
  axis.text.y = element_text(face="bold",size=15),
  axis.text.x = element_text(face="bold",size=15),
  panel.background = element_rect(fill="white"),
  axis.line = element_line(colour = "black", size = 1, linetype = "solid"))+
  scale_x_discrete(limit=blood_order) +
  annotate(geom="text", y=110, x=1.5, label=paste0("p=", round(pvalsBlood["USA300 dCrtN", "USA300"],4)), size=5)+
  annotate(geom="text", y=120, x=3.5, label=paste0("p=", round(pvalsBlood["SA925", "SA1088"],4)), size=5)

ggsave(BloodSurvivalPlot, file="Figures/Figure4/Figure4a_wholeblood.pdf", width=6, height=5)




# Figure 4b: whole blood survival of S. aureus 
##############################################

TTestNeutrophils = pairwise.t.test(neutrophils$survival, neutrophils$strain, pool.sd=F, p.adjust.method = "none")
pvalusNeutrophils= TTestNeutrophils$p.value

neutrophil_order<-c("USA300", "USA300 dCrtN","SA925","SA1088")
NeutrophilPlot = ggplot(neutrophils, aes(x=strain,y=survival))+geom_boxplot(fill="dodgerblue", width = 0.5,size=.25)+
  geom_jitter(width=.15)+
  labs(y="% Survival",x="strain")+
  theme(axis.title.x = element_text(face="bold",size=15),
        axis.title.y = element_text(face="bold",size=15),
        axis.text.y = element_text(face="bold",size=15),
        axis.text.x = element_text(face="bold",size=15),
        panel.background = element_rect(fill="white"),
        axis.line = element_line(colour = "black", 
                                 size = 1, linetype = "solid")) +
  scale_x_discrete(limit=neutrophil_order) + 
  annotate(geom="text",x=1.5, y=25, label=paste0("p=", round(pvalusNeutrophils["USA300 dCrtN", "USA300"], 4) ), size=5)+
  annotate(geom="text",x=3.5, y=28, label=paste0("p=", round(pvalusNeutrophils["SA925", "SA1088"], 4) ),size=5)

ggsave(NeutrophilPlot, file="Figures/Figure4/Figure4b_neutrophils.pdf", width=6, height=5)





# Figure 4(c): Neutrophil % in bl6
##################################

flowdata_bl6 =  flowdata %>% filter(MouseType=="C57BL6")
flowdata_bl6$condition = factor(flowdata_bl6$condition,levels=c("punch", "SA925", "SA1088" ))

# Same hypothesis testing Tori did originally 
Bl6nphil_punch = flowdata_bl6 %>% filter(condition=="punch")
Bl6nphil_SA1088 = flowdata_bl6 %>% filter(condition=="SA1088")
Bl6nphil_SA925 = flowdata_bl6 %>% filter(condition=="SA925")

NeutrophilsPunch1088 = t.test(Bl6nphil_punch$ly6g_percentage_cd11b,Bl6nphil_SA1088$ly6g_percentage_cd11b )
NeutrophilsPunch925 = t.test(Bl6nphil_punch$ly6g_percentage_cd11b,Bl6nphil_SA925$ly6g_percentage_cd11b )
NeutrophilsSA1088_925 = t.test(Bl6nphil_SA925$ly6g_percentage_cd11b,Bl6nphil_SA1088$ly6g_percentage_cd11b )

# Plot 
Bl6Neutrophils = ggplot(aes(x=condition, y=ly6g_percentage_cd11b), data=flowdata_bl6) + ylim(c(0,45))+
  geom_jitter(position=position_jitter(.1), size=2,color="dodgerblue")  + theme_classic() + 
  theme(axis.title.x = element_text(face="bold",size=15),axis.title.y = element_text(face="bold",size=15),
        axis.text.y = element_text(face="bold",size=15),axis.text.x = element_text(face="bold",size=15, angle=90, vjust=.5), 
        plot.title=element_text(face="bold", size=18,hjust=.5, vjust=1.5)) + 
  ylab("Ly6G+ Neutrophil Percentage (of CD11b cells)") + ggtitle("Neutrophil % of CD11b") +
  stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar", width = 0.3, color="black") + xlab("Condition") + 
  annotate(geom="text", y=42, x=1.5, label=paste0("p=",round(NeutrophilsPunch925$p.value, 4)),size=5) +
  annotate(geom="text", y=42, x=2.5, label=paste0("p=",round(NeutrophilsSA1088_925$p.value, 4)),size=5) + 
  annotate(geom="text", y=45, x=2, label=paste0("p=",round(NeutrophilsPunch1088$p.value, 4)),size=5)


ggsave(Bl6Neutrophils, file="Figures/Figure4/Figure4C_NeutrophilsBl6.pdf", height=6, width=5)


# Figure S4a
# Same comparison from Figure 4(a) in db/db mice
flowdata_dbdb =  flowdata %>% filter(MouseType=="dbdb" & condition!="none")
flowdata_dbdb$condition = factor(flowdata_dbdb$condition,levels=c("punch", "SA925", "SA1088" ))

# Same hypothesis testing Tori did originally 
dbdb_nphil_punch = flowdata_dbdb %>% filter(condition=="punch")
dbdb_nphil_SA1088 = flowdata_dbdb %>% filter(condition=="SA1088")
dbdb_nphil_SA925 = flowdata_dbdb %>% filter(condition=="SA925")
NeutrophilsDBDB_PunchSA1088 = t.test(dbdb_nphil_punch$ly6g_percentage_cd11b,dbdb_nphil_SA1088$ly6g_percentage_cd11b )
NeutrophilsDBDB_PunchSA925 = t.test(dbdb_nphil_punch$ly6g_percentage_cd11b,dbdb_nphil_SA925$ly6g_percentage_cd11b )
NeutrophilsDBDB_SA1088_SA925 = t.test(dbdb_nphil_SA925$ly6g_percentage_cd11b,dbdb_nphil_SA1088$ly6g_percentage_cd11b )

DBDBNeutrophils  = ggplot(aes(x=condition, y=ly6g_percentage_cd11b), data=flowdata_dbdb) +ylim(c(0,45))+
  geom_jitter(position=position_jitter(.1), size=2,color="dodgerblue")  + theme_classic() + 
  theme(axis.title.x = element_text(face="bold",size=15),axis.title.y = element_text(face="bold",size=15),
        axis.text.y = element_text(face="bold",size=15),axis.text.x = element_text(face="bold",size=15, angle=90, vjust=.5), 
        plot.title=element_text(face="bold", size=18,hjust=.5, vjust=1.5)) + 
  ylab("Ly6G+ Neutrophil Percentage (of CD11b cells)") + ggtitle("Neutrophil % of CD11b (db/db mice)") +
  stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar", width = 0.3, color="black") + xlab("Condition") + 
  annotate(geom="text", y=24, x=1.5, label=paste0("p=",round(NeutrophilsDBDB_PunchSA925$p.value, 4)),size=5) +
  annotate(geom="text", y=24, x=2.5, label=paste0("p=",round(NeutrophilsDBDB_SA1088_SA925$p.value, 4)),size=5) + 
  annotate(geom="text", y=26, x=2, label=paste0("p=",round(NeutrophilsDBDB_PunchSA1088$p.value, 4)),size=5)

ggsave(DBDBNeutrophils, file="Figures/Figure4/FigureS4A_NeutrophilsDBDB.pdf", width=5, height=6)



# Figure 4(d) : CD4+ T-cells % in bl6 
####################################

Bl6tcell_punch = flowdata_bl6 %>% filter(condition=="punch")
Bl6tcell_SA1088 = flowdata_bl6 %>% filter(condition=="SA1088")
Bl6tcell_SA925 = flowdata_bl6 %>% filter(condition=="SA925")

TCellsBl6Punch_1088 = t.test(Bl6tcell_punch$cd4_percentage_cd90.2,Bl6tcell_SA1088$cd4_percentage_cd90.2 )
TCellsBl6Punch_925 = t.test(Bl6tcell_punch$cd4_percentage_cd90.2,Bl6tcell_SA925$cd4_percentage_cd90.2 )
TCellsBl6_925_1088 = t.test(Bl6tcell_SA1088$cd4_percentage_cd90.2,Bl6tcell_SA925$cd4_percentage_cd90.2 )


Bl6CD4T = ggplot(aes(x=condition, y=cd4_percentage_cd90.2), data=flowdata_bl6) +ylim(c(0,40))+
  geom_jitter(position=position_jitter(.1), size=2,color="dodgerblue")  + theme_classic() + 
  theme(axis.title.x = element_text(face="bold",size=15),axis.title.y = element_text(face="bold",size=15),
        axis.text.y = element_text(face="bold",size=15),axis.text.x = element_text(face="bold",size=15, angle=90, vjust=.5), 
        plot.title=element_text(face="bold", size=18,hjust=.5, vjust=1.5)) + 
  ylab("CD4+ T cells Percentage (of CD90.2 cells)") + ggtitle("CD4+ % of CD90.2") +
  stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar", width = 0.3, color="black") + xlab("Condition") + 
  annotate(geom="text", y=35, x=1.5, label=paste0("p=",round(TCellsBl6Punch_925$p.value, 4)),size=5) +
  annotate(geom="text", y=35, x=2.5, label=paste0("p=",round(TCellsBl6_925_1088$p.value, 4)),size=5) + 
  annotate(geom="text", y=37, x=2, label=paste0("p=",round(TCellsBl6Punch_1088$p.value, 4)),size=5)

  

Bl6CD4Tyaxis = ggplot(aes(x=condition, y=cd4_percentage_cd90.2), data=flowdata_bl6) + ylim(c(0,45))+
  geom_jitter(position=position_jitter(.1), size=2,color="dodgerblue")  + theme_classic() + 
  theme(axis.title.x = element_text(face="bold",size=15),axis.title.y = element_text(face="bold",size=15),
        axis.text.y = element_text(face="bold",size=15),axis.text.x = element_text(face="bold",size=15, angle=90, vjust=.5), 
        plot.title=element_text(face="bold", size=18,hjust=.5, vjust=1.5)) + 
  ylab("CD4+ T cells Percentage (of CD90.2 cells)") + ggtitle("CD4+ % of CD90.2") +
  stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar", width = 0.3, color="black") + xlab("Condition")



ggsave(Bl6CD4T, file="Figures/Figure4/Figure4DTcellsBl6.pdf", width=5, height=6)


dbdb_tcell_punch = flowdata_dbdb %>% filter(condition=="punch")
dbdb_tcell_SA1088 = flowdata_dbdb %>% filter(condition=="SA1088")
dbdb_tcell_SA925 = flowdata_dbdb %>% filter(condition=="SA925")

cd4dbdbPunch_1088 = t.test(dbdb_tcell_punch$cd4_percentage_cd90.2,dbdb_tcell_SA1088$cd4_percentage_cd90.2 )
cd4dbdbSA925_1088 = t.test(dbdb_tcell_SA1088$cd4_percentage_cd90.2,dbdb_tcell_SA925$cd4_percentage_cd90.2 )
cd4dbdbPunch_925 = t.test(dbdb_tcell_SA925$cd4_percentage_cd90.2,dbdb_tcell_punch$cd4_percentage_cd90.2 )

DBDBcd4  = ggplot(aes(x=condition, y=cd4_percentage_cd90.2), data=flowdata_dbdb) + ylim(c(0,40))+
  geom_jitter(position=position_jitter(.1), size=2,color="dodgerblue")  + theme_classic() + 
  theme(axis.title.x = element_text(face="bold",size=15),axis.title.y = element_text(face="bold",size=15),
        axis.text.y = element_text(face="bold",size=15),axis.text.x = element_text(face="bold",size=15, angle=90, vjust=.5), 
        plot.title=element_text(face="bold", size=18,hjust=.5, vjust=1.5)) + 
  ylab("CD4+ T cells Percentage (of CD90.2 cells)") + ggtitle("CD4+ % of CD90.2 (db/db mice)") +
  stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar", width = 0.3, color="black") + xlab("Condition")+
  annotate(geom="text", y=35, x=1.5, label=paste0("p=",round(cd4dbdbPunch_925$p.value, 4)),size=5) +
  annotate(geom="text", y=35, x=2.5, label=paste0("p=",round(cd4dbdbSA925_1088$p.value, 4)),size=5) + 
  annotate(geom="text", y=37, x=2, label=paste0("p=",round(cd4dbdbPunch_1088$p.value, 4)),size=5)


ggsave(DBDBcd4, file="Figures/Figure4/FigureS4bTCellsdbdb.pdf", width=5, height=6)



# December 2022 Data
#####################

# crtN vs. LAC in db/db mice
############################
DBDB_Flow_CrtN = read.csv("Data/InVivoData/dbdb_Tori_December22.csv")
DBDB_Flow_CrtN = DBDB_Flow_CrtN %>% filter(condition!="none")

dbdb_neutrophil_punch = DBDB_Flow_CrtN %>% filter(condition=="punch")
dbdb_neutrophil_crtN = DBDB_Flow_CrtN %>% filter(condition=="crtn")
dbdb_neutrophil_LAC = DBDB_Flow_CrtN %>% filter(condition=="lac")

PunchCrtN = t.test(dbdb_neutrophil_punch$neutrophil_percent,dbdb_neutrophil_crtN$neutrophil_percent )
LAC_crtN = t.test(dbdb_neutrophil_crtN$neutrophil_percent,dbdb_neutrophil_LAC$neutrophil_percent )
Punch_LAC = t.test(dbdb_neutrophil_LAC$neutrophil_percent,dbdb_neutrophil_punch$neutrophil_percent  )



DBDB_CrtN  = ggplot(aes(x=condition, y=neutrophil_percent), data=DBDB_Flow_CrtN) + ylim(c(0,45))+
  geom_jitter(position=position_jitter(.1), size=2,color="dodgerblue")  + theme_classic() + 
  theme(axis.title.x = element_text(face="bold",size=15),axis.title.y = element_text(face="bold",size=15),
        axis.text.y = element_text(face="bold",size=15),axis.text.x = element_text(face="bold",size=15, angle=90, vjust=.5), 
        plot.title=element_text(face="bold", size=18,hjust=.5, vjust=1.5)) + 
  ylab("Neutrophils(%)") + ggtitle("Neutrophils (db/db mice)") +
  stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar", width = 0.3, color="black") + xlab("Condition") 
 
DBDB_CrtN$data$condition=factor(DBDB_CrtN$data$condition, levels=c("punch", "lac", "crtn"))
  
DBDB_CrtN = DBDB_CrtN + annotate(geom="text", y=41, x=1.5, label=paste0("p=",round(Punch_LAC$p.value, 4)),size=5) +
  annotate(geom="text", y=41, x=2.5, label=paste0("p=",round(LAC_crtN$p.value, 4)),size=5) + 
  annotate(geom="text", y=43, x=2, label=paste0("p=",round(PunchCrtN$p.value, 4)),size=5)




# crtN vs. LAC in bl6 mice 
############################

bl6_Flow = read.csv("Data/InVivoData/c57bl6_Tori_December22.csv")

bl6_neutrophil_punch = bl6_Flow %>% filter(condition=="punch")
bl6_neutrophil_crtN = bl6_Flow %>% filter(condition=="crtn")
bl6_neutrophil_LAC = bl6_Flow %>% filter(condition=="lac")



PunchCrtN_bl6 = t.test(bl6_neutrophil_punch$neutrophil_percent,bl6_neutrophil_crtN$neutrophil_percent )
LAC_crtN_bl6 = t.test(bl6_neutrophil_crtN$neutrophil_percent,bl6_neutrophil_LAC$neutrophil_percent )
Punch_LAC_bl6= t.test(bl6_neutrophil_LAC$neutrophil_percent,bl6_neutrophil_punch$neutrophil_percent )

bl6_Flow_LAC = bl6_Flow %>% filter(condition %in% c("punch", "crtn", "lac"))
bl6_CrtN  = ggplot(aes(x=condition, y=neutrophil_percent), data=bl6_Flow_LAC) + ylim(c(0,55))+
  geom_jitter(position=position_jitter(.1), size=2,color="dodgerblue")  + theme_classic() + 
  theme(axis.title.x = element_text(face="bold",size=15),axis.title.y = element_text(face="bold",size=15),
        axis.text.y = element_text(face="bold",size=15),axis.text.x = element_text(face="bold",size=15, angle=90, vjust=.5), 
        plot.title=element_text(face="bold", size=18,hjust=.5, vjust=1.5)) + 
  ylab("Neutrophils(%)") + ggtitle("Neutrophils (bl6 mice)") +
  stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar", width = 0.3, color="black") + xlab("Condition") 
bl6_CrtN$data$condition = factor(bl6_CrtN$data$condition, levels=c("punch", "lac", "crtn")) 
 
  
bl6_CrtN= bl6_CrtN+ annotate(geom="text", y=50, x=1.5, label=paste0("p=",round(Punch_LAC_bl6$p.value, 4)),size=5) +
  annotate(geom="text", y=50, x=2.5, label=paste0("p=",round(LAC_crtN_bl6$p.value, 4)),size=5) + 
  annotate(geom="text", y=54, x=2, label=paste0("p=",round(PunchCrtN_bl6$p.value, 4)),size=5)
ArrangedCrtComparison = grid.arrange(DBDB_CrtN,bl6_CrtN, ncol=2)
ggsave(ArrangedCrtComparison, file="Figures/UpToDateCrtN_LAC_Neutrophils.pdf", width=12, height=10)








# Pooling additional observations of SA1088 vs. SA925 in Bl6 (Tori says this is ok to do with cell type Pcts but not cell counts)

# Previous data:
##################

flowdata_bl6$neutrophil_percent = flowdata_bl6$ly6g_percentage_cd11b
Bl6nphil_punch = flowdata_bl6 %>% filter(condition=="punch") %>% select(condition, neutrophil_percent)
Bl6nphil_SA1088 = flowdata_bl6 %>% filter(condition=="SA1088")%>% select(condition, neutrophil_percent)
Bl6nphil_SA925 = flowdata_bl6 %>% filter(condition=="SA925")%>% select(condition, neutrophil_percent)

# Additional data points 
#########################
bl6_Flow = read.csv("Data/InVivoData/c57bl6_Tori_December22.csv")

bl6_neutrophil_punch = bl6_Flow %>% filter(condition=="punch") %>% select(condition, neutrophil_percent)
bl6_neutrophil_SA1088 = bl6_Flow %>% filter(condition=="1088") %>% select(condition, neutrophil_percent)
bl6_neutrophil_SA1088$condition="SA1088"
bl6_neutrophil_SA925 = bl6_Flow %>% filter(condition=="925") %>% select(condition, neutrophil_percent)
bl6_neutrophil_SA925$condition="SA925"



Combined_Bl6_punch = rbind(bl6_neutrophil_punch, Bl6nphil_punch)
Combined_Bl6_1088 = rbind(bl6_neutrophil_SA1088, Bl6nphil_SA1088)
Combined_Bl6_925 = rbind(bl6_neutrophil_SA925,Bl6nphil_SA925)

Black6_punch1088 = t.test(Combined_Bl6_1088$neutrophil_percent, Combined_Bl6_punch$neutrophil_percent)
Black6_punch925 = t.test(Combined_Bl6_925$neutrophil_percent, Combined_Bl6_punch$neutrophil_percent)
Black6_925_1088 =  t.test(Combined_Bl6_925$neutrophil_percent, Combined_Bl6_1088$neutrophil_percent)

fullDF_bl6 = rbind(Combined_Bl6_punch, Combined_Bl6_1088)
fullDF_bl6 = rbind(fullDF_bl6, Combined_Bl6_925)


bl6_neutrophils_clinical  = ggplot(aes(x=condition, y=neutrophil_percent), data=fullDF_bl6) + ylim(c(0,48))+
  geom_jitter(position=position_jitter(.1), size=2,color="dodgerblue")  + theme_classic() + 
  theme(axis.title.x = element_text(face="bold",size=15),axis.title.y = element_text(face="bold",size=15),
        axis.text.y = element_text(face="bold",size=15),axis.text.x = element_text(face="bold",size=15, angle=90, vjust=.5), 
        plot.title=element_text(face="bold", size=18,hjust=.5, vjust=1.5)) + 
  ylab("Neutrophils(%)") + ggtitle("Neutrophil % of CD11b (bl6 mice)") +
  stat_summary(fun = median, fun.min = median, fun.max = median, geom = "crossbar", width = 0.3, color="black") + xlab("Condition") 
bl6_neutrophils_clinical$data$condition = factor(bl6_neutrophils_clinical$data$condition, levels=c("punch", "SA925", "SA1088")) 

bl6_neutrophils_clinical = bl6_neutrophils_clinical + annotate(geom="text", y=43.5, x=1.5, label=paste0("p=",round(Black6_punch925$p.value, 4)),size=5) +
  annotate(geom="text", y=43.5, x=2.5, label=paste0("p=",round(Black6_925_1088$p.value, 4)),size=5) + 
  annotate(geom="text", y=46, x=2, label=paste0("p=",round(Black6_punch1088$p.value, 4)),size=5)

Clinical_Arranged = grid.arrange(bl6_neutrophils_clinical, DBDBNeutrophils,ncol=2)
ggsave(Clinical_Arranged, file="Figures/UpToDate_ClinicalIsolates_Flow.pdf", width=12, height=10)


# Corresponding CFU Counts
##########################
Black6CFUs = read.csv(file="Data/InVivoData/ear_cfu_count_c57BL6.csv")
dbdbCFUs = read.csv(file="Data/InVivoData/ear_cfu_count_dbdb_exp29.csv")

Black6CFUs = Black6CFUs %>% filter(!is.na(cfu_chromagar))
Black6CFUs = Black6CFUs %>% filter(condition %in% c("crtn", "lac", "1088", "925"))
dbdbCFUs = dbdbCFUs %>% filter(!is.na(cfu_chromagar))
dbdbCFUs = dbdbCFUs %>% filter(condition %in% c("crtn", "lac"))
bl6CFUplot = ggplot(Black6CFUs, aes(x=condition, y=log10(cfu_chromagar))) + geom_boxplot(fill="dodgerblue") + geom_jitter(height=0, width=.2) + theme_classic()+ stat_compare_means(method="t.test", comparisons = list(c("1088", "925"),c("lac", "crtn"))) +ylim(4,5)
dbdbCFUplot = ggplot(dbdbCFUs, aes(x=condition, y=log10(cfu_chromagar))) + geom_boxplot(fill="dodgerblue")+geom_jitter(height=0, width=.2) + theme_classic() + stat_compare_means(method="t.test",  comparisons = list(c("lac", "crtn")))+ylim(4,5)

compared = grid.arrange(bl6CFUplot, dbdbCFUplot, ncol=2,widths=2:1)
ggsave(compared,file="Figures/ToriCFUs.pdf", width=7, height=5)


Black6CFUs = read.csv(file="Data/InVivoData/ear_cfu_count_c57BL6.csv")
dbdbCFUs = read.csv(file="Data/InVivoData/ear_cfu_count_dbdb_exp29.csv")
Black6CFUs = Black6CFUs %>% filter(condition %in% c("crtn", "lac", "1088", "925"))
dbdbCFUs = dbdbCFUs %>% filter(condition %in% c("crtn", "lac"))
dbdbCFUs = dbdbCFUs %>% filter(!is.na(dilution_blood))
bl6totalCFUplot = ggplot(Black6CFUs, aes(x=condition, y=log10(cfu_ear_blood))) + geom_boxplot(fill="dodgerblue") + geom_jitter(height=0, width=.2) + theme_classic()+ stat_compare_means(method="t.test", comparisons = list(c("1088", "925"),c("lac", "crtn"))) +ylim(4,5)
dbdbCFUsPlot =  ggplot(dbdbCFUs, aes(x=condition, y=log10(cfu_ear_blood))) + geom_boxplot(fill="dodgerblue") + geom_jitter(height=0, width=.2) + theme_classic()+ stat_compare_means(method="t.test", comparisons = list(c("lac", "crtn"))) 
comparedTotal = grid.arrange(bl6totalCFUplot, dbdbCFUsPlot, ncol=2,widths=2:1)
ggsave(comparedTotal,file="Figures/ToriCFUsTotalBlood.pdf", width=7, height=5)
