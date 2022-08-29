# Figure 3
# adapted from ARM's scripts

woundsize= read.csv("Data/InVivoData/woundsizes.csv")

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


