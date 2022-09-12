
## Whole genome assembly & contamination checks


### Trim the raw reads of all 14 isolates in DFU141
```
# Trim the illumina reads
sh WholeGenomeSequencing/TrimIllumina.sh

# Trim the raw ONP reads 
sh WholeGenomeSequencing/RunPoreChopper.sh

```
### Run assembly with [unicycler](https://github.com/rrwick/Unicycler)


```
# Run unicycler in illumina mode (for SA929 + others' short-read based assemblies for contamination check)

# Run unicycler hybrid mode to perform hybrid assembly on all but SA929

```

### Contamination & coverage checks

```
# Run bowtie for trimmed illumina reads against Illumina-based assemblies 
Run_Bowtie_Assemblies.sh

# Use bowtie output to 
Coverage_Stats_Contigs.sh 

# .csv of  depths & lengths by contig
sh Merge_Coverage_Stats_DFU.sh



```
