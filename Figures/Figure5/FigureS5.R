library(ggplot2)
library(dplyr)


# RSBU KO PLOT
###############
RSBU = read.csv("Data/InVitroData/RsbU_stx_06.16.csv", header=T)

RSBU$Staphyloxanthin = 100*(RSBU$reading / 0.22566667)
RSBU = RSBU %>% mutate(Strain = if_else( name=="WT", "wildtype", paste0(paste0("\u0394", name) ))) 
RSBUPlot = ggplot(RSBU, aes(x=Strain, y=Staphyloxanthin)) + geom_boxplot(fill="darkorange", size=.25) + 
  geom_jitter(size=2, color="black") + theme_classic() + ylab("Staphyloxanthin Production(% USA300/LAC)") + ylim(0,110) 

pairwise.t.test(log(RSBU$Staphyloxanthin), RSBU$Strain, data= RSBU, p.adjust.method = "none")
ggsave(RSBUPlot, file="Figures/Figure5/FigureS5.pdf", width=6,height=7)
