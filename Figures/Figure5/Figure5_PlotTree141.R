# Amy Campbell
# 2022
# Plotting isolates from Patient 141

library("ape")
library("ggtree")
library("dplyr")
library("ggplot2")
library("ggtreeExtra")


TreeFilePath = "data/Patient141tree/RAxML_bestTree.RaxMLTree141Isolates"
UTDPhenotypes=read.csv("data/staphyloxanthin_paper_data.csv")
StaphIsolateDORNs$DORN=paste0("DORN", StaphIsolateDORNs$Doern.lab.bank.)
UTDPhenotypes = UTDPhenotypes %>% filter(patient==141)

TreeObject = ggtree::read.tree(TreeFilePath)
RootedTree = ape::root(TreeObject, "CC1_MW2",  resolve.root=T)
ggtree(RootedTree)

# Adjust Outgroup length
TreePlot <- ggtree::ggtree(RootedTree, size=.5, layout = "rectangular") + geom_tiplab(size=4, hjust=-.2)

# Scale down the outgroup by 1/10
TreePlot$data[TreePlot$data$label %in% c("CC1_MW2"), "x"] = (TreePlot$data[TreePlot$data$label %in% c("CC1_MW2"), "x"])*.05
TreePlot$data$DOERN = sapply(TreePlot$data$label, function(x) stringr::str_replace_all(string=x, pattern = "DORN", replacement = "SA"))
TreePlot + geom_treescale(label = "Mean # substitutions per site (across 2442 core genes)") 

UTDPhenotypes$DORN = paste0("DORN", UTDPhenotypes$DORN)
TreePlot$data$DORN = TreePlot$data$label
TreePlot$data = TreePlot$data %>% left_join(UTDPhenotypes, by="DORN")
TreePlot$data = TreePlot$data %>% mutate(VisitLabel = if_else( is.na(visit), "", paste0("Week ",visit*2 )
  
))
TreePlot$data$label = TreePlot$data$DOERN

TreePlot$data$VisitLabel = factor(TreePlot$data$VisitLabel, levels=c("Week 0","Week 2", "Week 4", "Week 6", "Week 8", "Week 10", "Week 12", "Week 16", "") )
UpdatedTreePlot = TreePlot + geom_tippoint(shape=15,size=5, aes(color=factor(VisitLabel))) + scale_color_manual(values=append(RColorBrewer::brewer.pal(n=8,name="Spectral"), "white"))  +
  geom_treescale(label = "Mean # substitutions per site (across 2442 core genes)") +
  guides(color=guide_legend(title="Collected")) + theme(legend.position="top") + xlim(0, .00002)
#  $label = sapply(TreePlot$data$visit, function(x) if_else(is.na(x), "", paste0( "(Week ", 2*x, ")") ))  
TreePlot$data = data.frame(TreePlot$data)
ggsave(UpdatedTreePlot, file="data/Patient141tree/Visualization141Isolates.pdf", width=11,height=8)  


UpdatedTreePlot + geom_fruit(aes(x=(staphyloxanthin)),fill="#B8860B", geom=geom_bar, stat="identity", orientation="y",offset=.3, pwidth=1)


# check for 

RoaryOutput =read.csv("data/Phylogeny2022Data/Patient141Roary/Roary/gene_presence_absence.csv")
RoaryOutput[,15:ncol(RoaryOutput)] <- if_else(RoaryOutput[,15:ncol(RoaryOutput)] =="", 0, 1)
RoaryOutput = RoaryOutput %>% select(-CC1_MW2)


RoaryOutput$rowsum = rowSums(RoaryOutput[,15:ncol(RoaryOutput)])

RoaryOutput %>% filter(DORN1000==0 & DORN1088==0 & DORN1038==0 & DORN1037==1 &  DORN1061==1 & DORN1194==1 & DORN880==1 & DORN881==1 & DORN882==1 & DORN925==1 & DORN929==1 & DORN933==1 & DORN976==1 & DORN999==1)
RoaryOutput %>% filter(DORN1000==1 & DORN1088==1 & DORN1038==1 & DORN1037==0 &  DORN1061==0 & DORN1194==0 & DORN880==0 & DORN881==0 & DORN882==0 & DORN925==0 & DORN929==0 & DORN933==0 & DORN976==0 & DORN999==0)

RoaryOutput %>% filter(DORN1088==1 & DORN925==0 )
RoaryOutput %>% filter(DORN925==1 & DORN1088==0 )


genomelistHigh = c(
"DORN1037",
"DORN1061",
"DORN1194",
"DORN880",
"DORN881",
"DORN882",
"DORN925",
"DORN929",
"DORN933",
"DORN976",
"DORN999")


genomelistHigh = c("DORN1000",
                   "DORN1037",
                   "DORN1038",
                   "DORN1061",
                   "DORN1088",
                   "DORN1194",
                   "DORN880",
                   "DORN881",
                   "DORN882",
                   "DORN925",
                   "DORN929",
                   "DORN933",
                   "DORN976",
                   "DORN999")

