
## Whole genome assembly & contamination checks


### Trim the raw reads of all 14 isolates in DFU141
```
# Trim the illumina reads
sh WholeGenomeSequencing/TrimIllumina.sh

# Trim the raw ONP reads 
sh WholeGenomeSequencing/RunPoreChopper.sh

```
### Run assembly with [unicycler](https://github.com/rrwick/Unicycler)

Run unicycler in illumina mode (for SA929 + others' short-read based assemblies for contamination check) as well as in hybrid mode for all but SA929 to get full assemblies

```
sh RunUnicycler.sh
```

### Contamination & coverage checks
Run bowtie for trimmed illumina reads against Illumina-based assemblies 
```
sh Run_Bowtie_Assemblies.sh
```
 Use bowtie output to summarize depths & breadths of reads by contig 
```
sh Coverage_Stats_Contigs.sh # calls coverage_stats.py
```
Make csv of  depths & lengths by contig (DFU_Covg_stats.tsv)

```
sh Merge_Coverage_Stats_DFU.sh
```

Plot coverage & sort contigs into two sets: one set that we’re letting be, and one that we’re blasting to check for contamination (outside of Median ± 1.5*IQR depth of coverage for the genome)
```
sh RunRScriptCovg.sh # calls CoverageAnalyses.R
sh SortContigs.sh # calls SortingContigsFate.py 
```

Run BLAST on the follow up contigs
```
sh BlastFollowupContigs.sh # uses BlastEnv and does NT blast on each ‘follow up’ contig
sh ReadBlastOutput.sh # Calls Add_SA_Contigs_Back.R to look at Blast output, see if any contigs blast to things other than S. Aureus 
```
Since no "follow-up" contigs blasted to non-S. aureus as a top hit, combine all >100 bp contigs back into a single file
```
sh CombineFilteredContigs.sh
```

Align reads against cleaned set of contigs, calculate coverage & call SNPs

```
TestAdmixture.sh # makes .bcf files from the alignments 
MakeSNPcsv.sh # calls getDP4Counts.py and outputs test_BCF_array.csv
```
Count SNP sites with >10X coverage ≥50 bp apart from one another and identify any genomes with more than 30 of these (based on [Raven et al. 2020](https://doi.org/10.1099/mgen.0.000354) )

