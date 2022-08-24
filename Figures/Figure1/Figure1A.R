# Amy Campbell
# Figure1 A Experimental setup
library("ggplot2")
library("gridExtra")
library("dplyr")
library("RColorBrewer")
library("forcats")

setwd("/Users/amycampbell/Desktop/GriceLabGit/Staphyloxanthin/")

FullData = read.csv("Data/InVitroData/staphyloxanthin_paper_updated.csv")

week_order <- c("2", "4", "6","8","10","12","14","16","Not Healed")

FullData$Healed = if_else(is.na(FullData$week_healed ), "Unhealed", "Healed")

FullData$week_collected<-(FullData$visit*2)

FullData = FullData %>% mutate(WeekHealed = if_else(is.na(week_healed),"Not Healed", as.character(week_healed)))


OrderByWeekHealed <- FullData %>% arrange(week_healed) 
OrderByWeekHealed$WeekHealed = factor(OrderByWeekHealed$WeekHealed)
OrderByWeekHealed$patient = factor(OrderByWeekHealed$patient)
IsolatesPlot = ggplot(data=OrderByWeekHealed, aes(x=week_collected, y =patient))+
  geom_tile(color="black", aes(fill=WeekHealed))+
  theme_classic()+ggtitle("S. aureus Isolates by Patient")+
  theme(plot.title = element_text(hjust = 0.5, face="bold", size=20))+
  labs(x="week collected",y="Patient")+
  theme(axis.title.x = element_text(face="bold",size=15))+
  theme(axis.title.y = element_text(face="bold",size=15))+
  theme(axis.text.y = element_blank())+
  theme(axis.text.x = element_text(face="bold",size=10))+  
  labs(fill="Week Healed") + scale_fill_manual(values=rev(RColorBrewer::brewer.pal(9, "Spectral"))) + scale_x_continuous(breaks=c(0,2,4,6,8,10,12,14,16,18,20,22,24,26)) 
IsolatesPlot$data$patient = factor(IsolatesPlot$data$patient, levels=unique(OrderByWeekHealed$patient))
IsolatesPlot$data$WeekHealed = factor(IsolatesPlot$data$WeekHealed, levels=week_order)

pdf(width=8.49, height= 6.555, file="Figures/Figure1/Fig1A.pdf")
IsolatesPlot
dev.off()

