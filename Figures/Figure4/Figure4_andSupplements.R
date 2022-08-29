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



