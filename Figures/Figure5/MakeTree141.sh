#!/bin/bash
# Amy Campbell
# Figure 5-- Making patient 141 isolate tree
# Run on cluster with -n 4 -M 10240 -R "rusage [mem=10240] span[hosts=1]"


# Run Prokka & roary
#####################
conda activate prokenv

# "genomes" folder in Data/ should contain just SA1088, SA1000, SA1038, SA1061,
# SA1037, SA999, SA976, SA933, SA929, SA882, SA881, SA880, CC1_MW2(outgroup)

gpath="../../Data/genomes/"
GFFs="../../Data/genomes/GFFfiles/"
justgffs="../../Data/genomes/GFFfiles/JustGFFs/"
roaryoutput="../../Data/genomes/Roary/"
trees="../../Data/Patient141tree/"
mkdir -p $trees 
mkdir -p $justgffs
mkdir -p $GFFs

extensionfasta="_Final.fasta"
gffext=".gff"
for genome in SA1000 SA1037 SA1038 SA1061 SA1088 SA1194 SA880 SA881 SA882 SA925 SA929 SA933 SA976 SA999 CC1_MW2; do 

   genomefile=$gpath$genome$extensionfasta
   prefixstring=$genome
   output=$GFFs$genome"/"

   prokka --outdir $output --force --prefix $prefixstring --genus Staphylococcus $genomefile
   echo $output$prefixstring$gffext
   echo $justgffs
   cp $output$prefixstring$gffext $justgffs

done

roary -e -p 4 -f $roaryoutput $justgffs*.gff


# Make the tree in RaxML
########################
conda activate TreeEnv

genealnstring="core_gene_alignment.aln"

coregenealignment=$roaryoutput$genealnstring

raxmlHPC-PTHREADS -T 4 -m GTRGAMMA -p 19104 -# 100 -s $coregenealignment -n RaxMLTree141Isolates
mv *RaxMLTree141Isolates* $trees 
