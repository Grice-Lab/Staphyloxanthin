#!/bin/bash
# Using Zhou et al.'s 2020 benchmarking study to look for genomes with >30 snps >50bp apart from one another 

source /home/acampbe/software/miniconda3/bin/activate BowtieEnv

Rscript AdmixtureOutput.R
