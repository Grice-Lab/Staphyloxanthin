# Figure 2 code for Staphyloxanthin paper
# Based on ARM and AEC's code 

library("ggplot2")
library("gridExtra")
library("dplyr")
library("RColorBrewer")
library("forcats")
library("ggpubr")


# Read in data
####################

# Replicate phenotype averages by patient, isolate, healing outcome
FullData  <-  read.csv("Data/InVitroData/staphyloxanthin_paper_updated.csv")

# Staphyloxanthin measures with replicates, 
patient141 <-read.csv("Data/InVitroData/patient_141_staphyloxanthin.csv")

# Staphyloxanthin measures for SA1088, SA925, USA300 LAC, delta-crtN-LAC 
StaphyloxanthinTestStrains =read.csv("Data/InVitroData/staphyloxanthin_by_strain.csv")

# Strain survival with hydrogen peroxide exposure, thymol inhibition 
StrainsH202 = read.csv("Data/InVitroData/hydrogenperoxide_survival_strain.csv")

# Strains staphyloxanthin thymol
StrainsThymolXanthin = read.csv("Data/InVitroData/staphyloxanthin_thymol.csv")

# Strains surival with polymixin 
StrainsPolymixin = read.csv("Data/InVitroData/polymixin_strain.csv")


# Figure 2A
###########
Staphyloxanthin_by_patient = ggplot(data=FullData, aes(x=factor(patient), y=staphyloxanthin, color=staphyloxanthin))+
  scale_color_gradient(low = "dodgerblue", high = "darkorange")+geom_point()+
 theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, face="bold", size=20), 
  axis.title.x = element_text(face="bold",size=15), 
  axis.title.y = element_text(face="bold",size=15),
  axis.text.y = element_text(face="bold",size=10),
  axis.text.x = element_text(face="bold",size=10, vjust=.5))+
  labs(color="Staphyloxanthin\nproduction(% 502A)", x="DFU",y="Staphyloxanthin")+
  theme(axis.text.x = element_text(angle = 90))

ggsave(Staphyloxanthin_by_patient, file="Figures/Figure2/Figure2A.pdf", width=9, height=5)

# Figure 2B
###########
IsolatesOrder = unique((patient141 %>% arrange(visit))$isolate)
patient141$weekcollect = factor(patient141$visit*2)

colorpalette = setdiff(rev(RColorBrewer::brewer.pal(10, "Spectral")), "#E6F598")
Patient141Plot <- ggplot(data=patient141, aes(x=isolate, y=staphyloxanthin)) + geom_boxplot(outlier.shape=NA,size=.1) + geom_jitter(width=.25, aes(color=weekcollect), size=2)+ scale_color_manual(values=colorpalette)+
  theme_classic()+
  theme(axis.title.x = element_text(face="bold",size=15), 
        axis.title.y = element_text(face="bold",size=15),
        axis.text.y = element_text(face="bold",size=15),
        axis.text.x = element_text(face="bold",size=10)) +
labs(color="Week collected", y="Staphyloxanthin production of biological replicate")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
Patient141Plot$data$isolate = factor(Patient141Plot$data$isolate, levels=IsolatesOrder)
ggsave(Patient141Plot, file="Figures/Figure2/Figure2B.pdf", width=8.175, height=5)


# Figure 2C
############
StaphyloxanthinTestStrains
strainorder <- c("USA300 LAC", "USA300 dCrtN", "SA925", "SA1088")

pairwise_Ts  = pairwise.t.test(StaphyloxanthinTestStrains$average, StaphyloxanthinTestStrains$strain, p.adjust.method = "none", pool.sd = F)


SA925_SA1088 = round((pairwise_Ts$p.value)["SA925", "SA1088"], 4)

SA925_LAC = round((pairwise_Ts$p.value)["USA300 LAC", "SA925"], 4)

LAC_CrtN =  round((pairwise_Ts$p.value)["USA300 LAC", "USA300 dCrtN"], 4)


StrainXanthinPlot = ggplot(StaphyloxanthinTestStrains, aes(x=strain,y=average))+geom_boxplot(fill="darkorange", width = 0.5, size=.25)+ 
  geom_jitter(width=.15)+
  labs(y="staphyloxanthin production (OD465 nm)",x="strain")+
  theme(plot.title = element_text(hjust = 0.5, face="bold", size=20),
        axis.title.x = element_text(face="bold",size=15),
        axis.title.y = element_text(face="bold",size=15),
        axis.text.y = element_text(face="bold",size=15),
        axis.text.x = element_text(face="bold",size=15),
        panel.background = element_rect(fill="white"),
        axis.line = element_line(colour = "black", size = 1, linetype = "solid")) +
  scale_x_discrete(limits = strainorder) +
  annotate(geom="text", label=paste0("p=", LAC_CrtN), y=.43, x=1.5, size=5)  + 
  annotate(geom="text", label=paste0("p=", SA925_LAC), y=.46, x=2, size=5) + 
  annotate(geom="text", label=paste0("p=", SA925_SA1088), y=.42, x=3.5, size=5)
        

ggsave(StrainXanthinPlot, width=7.5, height=6, file="Figures/Figure2/Figure2C.pdf")

# Figure 2D
###########
StrainsH202

peroxide_order<-c("USA300  LAC", "USA300 LAC dcrtN","SA925","SA1088","USA300  LAC +thymol","SA925 +thymol")
# For p values shown on ARM's figure 2D:
StrainsH202LACs = StrainsH202 %>% filter(strain %in% (c("USA300  LAC", "USA300 LAC dcrtN")))
StrainsH202Ttests = (pairwise.t.test(StrainsH202$sample, StrainsH202$strain, p.adjust.method = "none",pool.sd = F ))$p.value

LAC_crtN_h202 = round(StrainsH202Ttests["USA300 LAC dcrtN", "USA300  LAC"],4)
LAC_Thymol = round(StrainsH202Ttests["USA300  LAC +thymol", "USA300  LAC"],4)
LACcrtN_LACthymol = round(StrainsH202Ttests["USA300 LAC dcrtN","USA300  LAC +thymol"],4)
SA925_SA1088_h202 = round(StrainsH202Ttests["SA925","SA1088"],4)
SA925_thymol = round(StrainsH202Ttests["SA925 +thymol","SA925"],4)
SA1088_SA925thymol = round(StrainsH202Ttests["SA925 +thymol","SA1088"],4)

Strains_Survival_Thymol = ggplot(StrainsH202, aes(x=strain,y=sample))+geom_boxplot(fill="dodgerblue", width = 0.5, size=.25)+
  geom_jitter(width=.2)+
  labs(y="% Survival After Hydrogen Peroxide",x="strain")+
  theme(plot.title = element_text(hjust = 0.5, face="bold", size=20))+
  scale_x_discrete(limits = peroxide_order)+
  theme(axis.title.x = element_text(face="bold",size=15))+
  theme(axis.title.y = element_text(face="bold",size=15))+
  theme(axis.text.y = element_text(face="bold",size=10))+
  theme(axis.text.x = element_text(face="bold",size=8))+
  theme(panel.background = element_rect(fill="white"))+
  theme( axis.line = element_line(colour = "black", 
                                  size = 1, linetype = "solid"))+
  annotate(geom="text", y=138, x=1.5, label=paste0("p=", LAC_crtN_h202), size=5) +
  annotate(geom="text", y=145, x=3, label=paste0("p=", LAC_Thymol), size=5) +
  annotate(geom="text", y=138, x=3, label=paste0("p=", LACcrtN_LACthymol), size=5) + 
  annotate(geom="text", y=131, x=3.5, label=paste0("p=", SA925_SA1088_h202), size=5)+
  annotate(geom="text", y=124, x=4.5, label=paste0("p=", SA925_thymol), size=5)+
  annotate(geom="text", y=117, x=5, label=paste0("p=", SA1088_SA925thymol), size=5)

ggsave(Strains_Survival_Thymol, file="Figures/Figure2/Figure2D.pdf", width=8, height=6)


# Figure 2E
############

thymol_order <- c("USA300 LAC", "USA300 dCrtN", "SA925", "SA1088")

USA300fig2e = StrainsThymolXanthin %>% filter(strain=="USA300 LAC")
delta_crtNfig2e = StrainsThymolXanthin %>% filter(strain=="USA300 dCrtN")
SA925fig2e = StrainsThymolXanthin %>% filter(strain=="SA925")
SA1088fig2e = StrainsThymolXanthin %>% filter(strain=="SA1088")

t_test_USA300fig2e = t.test(USA300fig2e$staphyloxanthin ~ USA300fig2e$thymol)
t_test_delta_crtNfig2e = t.test(delta_crtNfig2e$staphyloxanthin ~ delta_crtNfig2e$thymol)
t_test_1088_fig2e = t.test(SA1088fig2e$staphyloxanthin ~ SA1088fig2e$thymol)
t_test_925_fig2e = t.test(SA925fig2e$staphyloxanthin ~ SA925fig2e$thymol)

pvalUSA300fig2e = round(t_test_USA300fig2e$p.value,4)
pvalDeltaCrtNfig2e = round( t_test_delta_crtNfig2e$p.value,4)
pval1088fig2e = round(t_test_1088_fig2e$p.value, 4)
pval925fig2e = "1.607e-05"

ThymolStrains = ggplot(StrainsThymolXanthin, aes(x=strain, y=staphyloxanthin, fill=thymol))+
  geom_jitter( position=position_jitterdodge(), size=1)+
  geom_boxplot( size=.2)+
  scale_fill_manual(values = c("darkorange", 'magenta')) +
  labs(x="strain", y="staphyloxanthin production (OD465 nm)") +
  theme_classic()+
  theme(axis.title.x = element_text(face="bold",size=15),
        axis.title.y = element_text(face="bold",size=15),
        axis.text.y = element_text(face="bold",size=15),
        axis.text.x = element_text(face="bold",size=15))+
  scale_x_discrete(limits = thymol_order) + 
  annotate(geom="text", label=paste0("p=", pvalUSA300fig2e), y=.45, x=1, size=5)+
  annotate(geom="text", label=paste0("p=", pvalDeltaCrtNfig2e), y=.45, x=2, size=5)+
  annotate(geom="text", label=paste0("p=", pval925fig2e), y=.45, x=3, size=5) + 
  annotate(geom="text", label=paste0("p=", pval1088fig2e), y=.45,x=4, size=5)

ggsave(ThymolStrains, width=8, height=6, file="Figures/Figure2/Figure2e.pdf")

      


# Figure 2F
############

USA300StrainsComparePoly = StrainsPolymixin %>% filter(strain %in% c("USA300", "USA300 dCrtN"))
ClinicalStrainsComparePoly = StrainsPolymixin %>% filter(strain %in% c("SA925", "SA1088"))

polymixin_t_USA300 = t.test(USA300StrainsComparePoly$survival ~ USA300StrainsComparePoly$strain)
polymixin_t_clinical = t.test(ClinicalStrainsComparePoly$survival ~ ClinicalStrainsComparePoly$strain)
p_polymixin300 = round(polymixin_t_USA3ThymolStrains00$p.value, 4)
p_polymixinClinical = round(polymixin_t_clinical$p.value, 4)

PolymixinPlot = ggplot(StrainsPolymixin, aes(x=strain,y=survival))+geom_boxplot(fill="dodgerblue", width = 0.5, size=.2 )+
  geom_jitter(width=.15)+ theme_classic() +
  labs(y="% Survival With Polymixin" ,x="strain")+
  theme(axis.title.x = element_text(face="bold",size=15),
        axis.title.y = element_text(face="bold",size=15),
        axis.text.y = element_text(face="bold",size=15),
        axis.text.x = element_text(face="bold",size=15),
        panel.background = element_rect(fill="white"), 
        axis.line = element_line(colour = "black", 
                                 size = 1, linetype = "solid"))+
  scale_x_discrete(limit=c("USA300", "USA300 dCrtN","SA925","SA1088")) +
  annotate(geom="text", x=1.5, y=125, label=paste0("p=", p_polymixin300), size=5) + 
  annotate(geom="text", x=3.5, y=125, label=paste0("p=", p_polymixinClinical),size=5)

ggsave(PolymixinPlot, width=7, height=6,file="Figures/Figure2/figure2F.pdf")

# Figure S2
###########
Palette5 = (RColorBrewer::brewer.pal(10, "Spectral"))[c(4, 1, 9, 8, 10)]

Palette5 = (RColorBrewer::brewer.pal(10, "Spectral"))[c(4, 1, 9, 8)]
Palette5 = (RColorBrewer::brewer.pal(10, "Spectral"))[c( 1, 9, 8, 10)]
#growth curves
growth<-read.csv("Data/InVitroData/growth_curves.csv")
growth
summary(growth)
p <- ggplot(data = growth, aes(x = time, y = average, group = strain))
p+geom_point()
growthcurve = p + geom_line(aes(color=strain),size=1)+
  geom_point(aes(color=strain),size=5)+
  scale_color_manual(values=Palette5)+
  theme_classic()+ggtitle("Growth Curve")+
  theme(plot.title = element_text(hjust = 0.5, face="bold", size=20))+
  labs(x="time (hours)",y="OD600")+
  theme(axis.title.x = element_text(face="bold",size=15))+
  theme(axis.title.y = element_text(face="bold",size=15))+
  theme(axis.text.y = element_text(face="bold",size=15))+
  theme(axis.text.x = element_text(face="bold",size=15))+
  geom_errorbar(aes(ymin =average-stdev, ymax =average+stdev),
                width = 0.2, position = position_dodge(0.5))
ggsave(growthcurve, file="Figures/Figure2/FigureS2.pdf", width=8,height=7)
