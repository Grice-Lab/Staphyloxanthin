

# Figure 3B
# only males were measured at Day 21 -- was this because of
#wounds now


wounds = read.csv("Data/InVivoData//woundsizes.csv")




wound_order <- c("PBS","LAC", "LAC_dCrtN_1", "SA925", "SA1088")

Woundplot = ggplot(data=wounds, aes(x=group, y=day_14_change))+geom_boxplot(size=.25, outlier.shape = NA)+
  geom_jitter(aes(color=group, shape=sex),size=2, width=.2)+
  labs(x="group", y="% original wound size at day 14")+ theme_classic()+
  theme(axis.title.x = element_text(face="bold",size=15))+
  theme(axis.title.y = element_text(face="bold",size=15))+
  theme(axis.text.y = element_text(face="bold",size=15))+
  theme(axis.text.x = element_text(face="bold",size=15))+
  theme(plot.title = element_text(hjust = 0.5, face="bold", size=20))+
  scale_color_manual(values=c("dodgerblue","cyan4","magenta","darkorange","blueviolet"))+
  scale_x_discrete(limits = wound_order)
ggsave(Woundplot, file="Figures/Figure3/Figure3b.pdf", width=9.6, height=4.88)

WilcoxResults = pairwise.wilcox.test(wounds$day_14_change, wounds$group, p.adjust.method = "BH", pool.sd=F)
WilcoxResultsPValues = WilcoxResults$p.value


pvalPBS_LAC = "6.165e-06"
pvalPBS_LAC = "6.165e-06"

round(WilcoxResultsPValues["PBS", "LAC"])
