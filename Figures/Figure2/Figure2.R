# Figure 2 code for Staphyloxanthin paper
# Based on ARM and AEC's code 

library("ggplot2")
library("gridExtra")
library("dplyr")
library("RColorBrewer")
library("forcats")
library("ggpubr")

WeekPalette =RColorBrewer::brewer.pal(11, "Spectral")

darkpurple =RColorBrewer::brewer.pal(11, "PuOr")[10]
WeekPalette = c("#8B0000", WeekPalette)

WeekPalette = c("black", WeekPalette)
WeekPalette = append(WeekPalette, darkpurple)
colormap=data.frame(colors=WeekPalette, week=c(0,2,4,6,8,10,12,14,16,18,20,22,24,26))

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
FullData$week = 2*(FullData$visit)
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

Staphyloxanthin_by_patient = ggplot(data=FullData, aes(x=factor(patient), y=staphyloxanthin, color=factor(week)))+
  scale_color_manual(values = colormap$colors)+geom_point()+
  theme_classic()+
  theme(plot.title = element_text(hjust = 0.5, face="bold", size=20), 
        axis.title.x = element_text(face="bold",size=15), 
        axis.title.y = element_text(face="bold",size=15),
        axis.text.y = element_text(face="bold",size=10),
        axis.text.x = element_text(face="bold",size=10, vjust=.5))+
  labs(color="Week collected", x="DFU",y="Staphyloxanthin")+
  theme(axis.text.x = element_text(angle = 90)) + geom_vline(xintercept = seq(.5, 61,1), linewidth=.05)

ggsave(Staphyloxanthin_by_patient, file="Figures/Figure2/Figure2A_2023.pdf", width=9, height=3)

# Figure 2B
###########
IsolatesOrder = unique((patient141 %>% arrange(visit))$isolate)
patient141$week = factor(patient141$visit*2)
colormap$week = factor(colormap$week)
#colorpalette = setdiff(rev(RColorBrewer::brewer.pal(10, "Spectral")), "#E6F598")
colorsdf =patient141 %>% select(week) %>% unique() %>% left_join(colormap, by="week")
Patient141Plot <- ggplot(data=patient141, aes(x=isolate, y=staphyloxanthin)) + geom_boxplot(outlier.shape=NA,size=.1) + geom_jitter(width=.25, aes(color=week), size=2,height=0)+ scale_color_manual(values=colorsdf$colors)+
  theme_classic()+
  theme(axis.title.x = element_text(face="bold",size=15), 
        axis.title.y = element_text(face="bold",size=15),
        axis.text.y = element_text(face="bold",size=15),
        axis.text.x = element_text(face="bold",size=10)) +
labs(color="Week collected", y="Staphyloxanthin production of biological replicate")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
Patient141Plot$data$isolate = factor(Patient141Plot$data$isolate, levels=IsolatesOrder)
ggsave(Patient141Plot, file="Figures/Figure2/Figure2B_2023.pdf", width=7, height=5)


# Figure 2C
############
StaphyloxanthinTestStrains
strainorder <- c("USA300 LAC", "USA300 dCrtN", "SA925", "SA1088")

pairwise_Ts  = pairwise.t.test(StaphyloxanthinTestStrains$average, StaphyloxanthinTestStrains$strain, p.adjust.method = "none", pool.sd = F)


SA925_SA1088 = round((pairwise_Ts$p.value)["SA925", "SA1088"], 4)

SA925_LAC = round((pairwise_Ts$p.value)["USA300 LAC", "SA925"], 4)

LAC_CrtN =  round((pairwise_Ts$p.value)["USA300 LAC", "USA300 dCrtN"], 4)


StrainXanthinPlot = ggplot(StaphyloxanthinTestStrains, aes(x=strain,y=average))+geom_boxplot(fill="darkorange", width = 0.5, size=.25)+ 
  geom_jitter(width=.15,height=0)+
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
        

ggsave(StrainXanthinPlot, width=7.5, height=3.5, file="Figures/Figure2/Figure2C_2023.pdf")

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
  geom_jitter(width=.2,height=0)+
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

ggsave(Strains_Survival_Thymol, file="Figures/Figure2/Figure2D_2023.pdf", width=8, height=3.5)


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

ggsave(ThymolStrains, width=9, height=3.5, file="Figures/Figure2/Figure2e_2023.pdf")

      


# Figure 2F
############

USA300StrainsComparePoly = StrainsPolymixin %>% filter(strain %in% c("USA300", "USA300 dCrtN"))
ClinicalStrainsComparePoly = StrainsPolymixin %>% filter(strain %in% c("SA925", "SA1088"))

polymixin_t_USA300 = t.test(USA300StrainsComparePoly$survival ~ USA300StrainsComparePoly$strain)
polymixin_t_clinical = t.test(ClinicalStrainsComparePoly$survival ~ ClinicalStrainsComparePoly$strain)

p_polymixin300 = round(polymixin_t_USA300$p.value, 4)
p_polymixinClinical = round(polymixin_t_clinical$p.value, 4)

PolymixinPlot = ggplot(StrainsPolymixin, aes(x=strain,y=survival))+geom_boxplot(fill="dodgerblue", width = 0.5, size=.2 )+
  geom_jitter(width=.15,height=0)+ theme_classic() +
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

ggsave(PolymixinPlot, width=8, height=4,file="Figures/Figure2/figure2F.pdf")


# Figure 2H
############
# staphyloxanthin production by complemented transposon mutant 

CompSTX = read.csv("Data/InVitroData/Complement_STX.csv")

CompSTX_je2_crtn = CompSTX %>% filter(Strain %in% c("JE2", "JE2_deltaCrtN"))
JE2crtNt = t.test(OD450_600 ~ Strain,data=CompSTX_je2_crtn)
JE2crtN_p = round(JE2crtNt$p.value,4)

CompSTX_crtn_empty = CompSTX %>% filter(Strain %in% c( "JE2_deltaCrtN", "JE2_deltaCrtN_pepSA5"))
crtN_empty_T = t.test(OD450_600 ~ Strain,data=CompSTX_crtn_empty)
crtN_empty_p = round(crtN_empty_T$p.value,4)

CompSTX_JE2_Complement = CompSTX %>% filter(Strain %in% c( "JE2", "JE2_deltaCrtN_pepSA5-CrtMN"))
complement_JE2_T = t.test(OD450_600 ~ Strain,data=CompSTX_JE2_Complement)
complement_JE2_p = round(complement_JE2_T$p.value, 4)

CompSTX_empty_complement = CompSTX %>% filter(Strain %in% c( "JE2_deltaCrtN_pepSA5", "JE2_deltaCrtN_pepSA5-CrtMN"))
complement_empty_T = t.test(OD450_600 ~ Strain,data=CompSTX_empty_complement)
complement_empty_p = round(complement_empty_T$p.value, 4)


compstx_plot = ggplot(CompSTX, aes(x=Strain, y=OD450_600)) + 
  geom_boxplot(fill="darkorange") + theme_classic() + geom_jitter(height=0)+
  annotate(geom="text", x=1.5, y=1.6, label=paste0("p=",JE2crtN_p ), size=5) + 
  annotate(geom="text", x=2.5, y=.5, label=paste0("p=",crtN_empty_p ), size=5) + 
  annotate(geom="text", x=2.5, y=1.6, label=paste0("p=",complement_JE2_p ), size=5) + 
  annotate(geom="text", x=3.5, y=1.7, label=paste0("p=",complement_empty_p ), size=5) 
ggsave(compstx_plot, width=8, height=4,file="Figures/Figure2/figure2H.pdf")


# Figure 2I
###########
# 5-41 = JE2
# 7-08 = ∆crtN-pEPS5a
# 7-10 = ∆crtN-pEP-CrtMN

H2O2 =read.csv2('Data/InVitroData/H2O2_Complement_CFUs.txt', sep="\t")
H2O2$RelDPBS = sapply(H2O2$RelDPBS, function(x) as.numeric(as.character(x)))


# choose middle dilution for .30% concentration
H2O2 = H2O2 %>% filter( (H2O2_Per=="0.30%")) %>% filter(DilutionID=="D5")
shapiro.test(H2O2$RelDPBS)


H2O2 = H2O2 %>% mutate(Strain = case_when(SaStrain=="EGM5-41" ~ "JE2",
                                          SaStrain=="EGM7-08" ~ "crtNpEPSA5", 
                                          SaStrain=="EGM7-10"~"crtNpEP-CrtMN"))

# So that we can pair by experiment 
H2O2 = H2O2 %>% arrange(SaStrain,Experiment)

H2O2$PctSurvival = H2O2$RelDPBS*100
CompareJE2Empty = H2O2 %>% filter(Strain %in% c("JE2", "crtNpEPSA5"))
CompareJE2Empty_p = round((t.test(PctSurvival~Strain,data=CompareJE2Empty, paired=T))$p.value, 4)

CompareEmpty_Complement = H2O2 %>% filter(Strain %in% c("crtNpEPSA5", "crtNpEP-CrtMN"))
CompareEmpty_Complement_p = round((t.test(PctSurvival~Strain,data=CompareEmpty_Complement, paired=T))$p.value, 4)

CompareJE2_Complement=H2O2 %>% filter(Strain %in% c("JE2", "crtNpEP-CrtMN"))
CompareJE2_Complement_p = round((t.test(PctSurvival~Strain,data=CompareJE2_Complement, paired=T))$p.value, 4)

h2o2Complement=ggplot(H2O2, aes(x=factor(Strain), y=PctSurvival)) + geom_boxplot(fill="dodgerblue") + geom_jitter(height=0,width=.2) + theme_classic() + scale_x_discrete(limits=c("JE2","crtNpEPSA5", "crtNpEP-CrtMN"))
h2o2Complement = h2o2Complement + annotate(geom="text", x=1.5, y=150, label=paste0("p=",CompareJE2Empty_p ), size=5) + 
  annotate(geom="text", x=2.5, y=200, label=paste0("p=",CompareEmpty_Complement_p ), size=5) + 
  annotate(geom="text", x=2, y=210, label=paste0("p=",CompareJE2_Complement_p ), size=5) 
ggsave(h2o2Complement, width=8, height=4,file="Figures/Figure2/figure2I.pdf")

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
