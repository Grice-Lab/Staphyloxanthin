#!/bin/bash
# Align and compare whole genomes from patient 141
# specifically in this case compare each genome to 925

# NOTE: change this to wherever you've cloned the repo. You must have downloaded the NCBI genomes to Staphyloxanthin/Data/AllGenomes/

WD=Staphyoxanthin

source /home/acampbe/software/miniconda3/bin/activate Patient141

mkdir -p $WD/Data/WGSData/NucDiffOutput/


nucdiff $WD/Data/AllGenomes/SA925._Final.fasta \
        $WD/Data/AllGenomes/SA1000._Final.fasta \
        $WD/Data/WGSData/NucDiffOutput/NucDiff_925_1000 SA925_SA1000Comparison

nucdiff $WD/Data/AllGenomes/SA925._Final.fasta \
        $WD/Data/AllGenomes/SA976._Final.fasta \
        $WD/Data/WGSData/NucDiffOutput/NucDiff_925_976 SA925_SA976Comparison

nucdiff $WD/Data/AllGenomes/SA925._Final.fasta \
        $WD/Data/AllGenomes/SA999._Final.fasta \
        $WD/Data/WGSData/NucDiffOutput/NucDiff_925_999 SA925_SA999Comparison

nucdiff $WD/Data/AllGenomes/SA925._Final.fasta \
        $WD/Data/AllGenomes/SA1194._Final.fasta \
        $WD/Data/WGSData/NucDiffOutput/NucDiff_925_1194 SA925_SA1194Comparison

nucdiff $WD/Data/AllGenomes/SA925._Final.fasta \
        $WD/Data/AllGenomes/SA1061._Final.fasta \
        $WD/Data/WGSData/NucDiffOutput/NucDiff_925_1061 SA925_SA1061Comparison

nucdiff $WD/Data/AllGenomes/SA925._Final.fasta \
        $WD/Data/AllGenomes/SA1088._Final.fasta \
        $WD/Data/WGSData/NucDiffOutput/NucDiff_925_1088 SA925_SA1088Comparison

nucdiff $WD/Data/AllGenomes/SA925._Final.fasta \
        $WD/Data/AllGenomes/SA1038._Final.fasta \
        $WD/Data/WGSData/NucDiffOutput/NucDiff_925_1038 SA925_SA1038Comparison

nucdiff $WD/Data/AllGenomes/SA925._Final.fasta \
        $WD/Data/AllGenomes/SA881._Final.fasta \
        $WD/Data/WGSData/NucDiffOutput/NucDiff_925_881 SA925_SA881Comparison

nucdiff $WD/Data/AllGenomes/SA925._Final.fasta \
        $WD/Data/AllGenomes/SA880._Final.fasta \
        $WD/Data/WGSData/NucDiffOutput/NucDiff_925_880 SA925_SA880Comparison

nucdiff $WD/Data/AllGenomes/SA925._Final.fasta \
        $WD/Data/AllGenomes/SA1037._Final.fasta \
        $WD/Data/WGSData/NucDiffOutput/NucDiff_925_1037 SA925_SA1037Comparison


nucdiff $WD/Data/AllGenomes/SA925._Final.fasta \
        $WD/Data/AllGenomes/SA933._Final.fasta \
        $WD/Data/WGSData/NucDiffOutput/NucDiff_925_933 SA925_SA933Comparison


nucdiff $WD/Data/AllGenomes/SA925._Final.fasta \
        $WD/Data/AllGenomes/SA882._Final.fasta \
        $WD/Data/WGSData/NucDiffOutput/NucDiff_925_882 SA925_SA882Comparison


