#!/bin/bash
# Amy Campbell
# Figure 5-- Making patient 141 isolate tree
# Run on cluster with -n 4 -M 10240 -R "rusage [mem=10240] span[hosts=1]"


# Run Prokka & roary
#####################
conda activate prokenv

gpath="../../Data/genomes/"
GFFs="../../Data/genomes/GFFfiles/"
justgffs="../../Data/genomes/GFFfiles/JustGFFs/"
roaryoutput="../../Data/genomes/Roary/"

mkdir -p $justgffs
mkdir -p $GFFs

extensionfasta="_Final.fasta"
gffext=".gff"
for genome in DORN1000 DORN1037 DORN1038 DORN1061 DORN1088 DORN1194 DORN880 DORN881 DORN882 DORN925 DORN929 DORN933 DORN976 DORN999 CC1_MW2; do 

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
