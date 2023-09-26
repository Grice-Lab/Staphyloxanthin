#!/bin/bash
# Align and compare whole genomes from patient 141
# specifically in this case compare each genome to 925

# NOTE: change this to wherever you've cloned the repo. You must have downloaded the NCBI genomes to Staphyloxanthin/Data/AllGenomes/

WD=Staphyoxanthin/

source /home/acampbe/software/miniconda3/bin/activate Patient141

mkdir -p /home/acampbe/DFU/data/WGS_2020/ONP/NucDiff925/


nucdiff /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/SA925._Final.fasta \
        /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/SA1000._Final.fasta \
        /home/acampbe/DFU/data/WGS_2020/ONP/NucDiff925/NucDiff_925_1000 SA925_SA1000Comparison

nucdiff /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/SA925._Final.fasta \
        /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/SA976._Final.fasta \
        /home/acampbe/DFU/data/WGS_2020/ONP/NucDiff925/NucDiff_925_976 SA925_SA976Comparison

nucdiff /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/SA925._Final.fasta \
        /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/SA999._Final.fasta \
        /home/acampbe/DFU/data/WGS_2020/ONP/NucDiff925/NucDiff_925_999 SA925_SA999Comparison

nucdiff /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/SA925._Final.fasta \
        /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/SA1194._Final.fasta \
        /home/acampbe/DFU/data/WGS_2020/ONP/NucDiff925/NucDiff_925_1194 SA925_SA1194Comparison

nucdiff /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/SA925._Final.fasta \
        /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/SA1061._Final.fasta \
        /home/acampbe/DFU/data/WGS_2020/ONP/NucDiff925/NucDiff_925_1061 SA925_SA1061Comparison

nucdiff /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/SA925._Final.fasta \
        /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/SA1088._Final.fasta \
        /home/acampbe/DFU/data/WGS_2020/ONP/NucDiff925/NucDiff_925_1088 SA925_SA1088Comparison

nucdiff /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/SA925._Final.fasta \
        /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/SA1038._Final.fasta \
        /home/acampbe/DFU/data/WGS_2020/ONP/NucDiff925/NucDiff_925_1038 SA925_SA1038Comparison

nucdiff /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/SA925._Final.fasta \
        /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/SA881._Final.fasta \
        /home/acampbe/DFU/data/WGS_2020/ONP/NucDiff925/NucDiff_925_881 SA925_SA881Comparison

nucdiff /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/SA925._Final.fasta \
        /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/SA880._Final.fasta \
        /home/acampbe/DFU/data/WGS_2020/ONP/NucDiff925/NucDiff_925_880 SA925_SA880Comparison

nucdiff /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/SA925._Final.fasta \
        /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/SA1037._Final.fasta \
        /home/acampbe/DFU/data/WGS_2020/ONP/NucDiff925/NucDiff_925_1037 SA925_SA1037Comparison


nucdiff /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/SA925._Final.fasta \
        /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/SA933._Final.fasta \
        /home/acampbe/DFU/data/WGS_2020/ONP/NucDiff925/NucDiff_925_933 SA925_SA933Comparison


nucdiff /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/SA925._Final.fasta \
        /home/acampbe/DFU/data/WGS_2020/FinalSetStaphIsolates/SA882._Final.fasta \
        /home/acampbe/DFU/data/WGS_2020/ONP/NucDiff925/NucDiff_925_882 SA925_SA882Comparison


