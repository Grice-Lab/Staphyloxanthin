# Amy  Campbell
# Updated 01/2022 -- adapted for 09/22 staphyloxanthin paper
# Sorting contigs for BLAST followup or keep (based on coverage & length as analyzed in ContigDepthAnalyses.R)

# Iterate through all the .fasta files in the cleaned assemblies
# For each one,
# Find the subset of lines in ContigFilteringFates.csv with Genome == (that file's basename but with _contigs.fasta removed)
# With seqio, open the fasta as a dictionary

# Rename all the records so that they are just the first (if that's not how they're keyed anyway?)
# Name the file with the "keep" contigs <GENOMENAME>_CleanContigs.fasta
# Name the file with the "follow up" contigs <GENOMENAME>_FollowUpContigs.fasta

#
from Bio import SeqIO
import pandas
import os
import sys

# arg 1 should be path to assemblies
# arg 2 should be path to contig fates
# arg 3 should be path to output cleaned contigs
# arg 4 should be path to output follow up contigs
assemblypath = str(sys.argv[1])
contigfatespath = str(sys.argv[2])
cleancontigspath=str(sys.argv[3])
follupcontigspath=str(sys.argv[4])

contigfates = pandas.read_csv(contigfatespath)


# Read in each fasta file (uncleaned contigs)
##############################################
for genomefile in os.listdir(assemblypath):
	
    genomename = genomefile.replace("_contigs.fasta", "")

    print("Now processing " + genomename + "...")

    # Read in this fasta file as a SeqIO dictionary where keys are contigID and values are the remaining info + sequence
    filepath = assemblypath + str(genomefile)
    SeqDictionary = SeqIO.to_dict(SeqIO.parse(filepath, "fasta"))

    # Subset the contigfates dataframe to include only info about the genome we're dealing with
    reduced_Fates = contigfates[contigfates["Genome"] == str(genomename)]

    reduced_Fates_followup = reduced_Fates[reduced_Fates["SortContig"]=="FollowUp"]
    reduced_Fates_keep = reduced_Fates[reduced_Fates["SortContig"]=="keep"]

    key_list_followup = list(map(str, list(reduced_Fates_followup["ContigID"].values)))
    key_list_keep = list(map(str, list(reduced_Fates_keep["ContigID"].values)))

    SeqDictionary_followup = {k: SeqDictionary[k] for k in key_list_followup}
    SeqDictionary_keep = {k: SeqDictionary[k] for k in key_list_keep}

    followup_out = open(str(follupcontigspath +genomename + "_FollowUpContigs.fasta"), "w")
    SeqIO.write(SeqDictionary_followup.values(), followup_out, "fasta")
    followup_out.close()

    keep_out = open(str(cleancontigspath  + genomename + "_CleanContigs.fasta"), "w")
    SeqIO.write(SeqDictionary_keep.values(), keep_out, "fasta")
    keep_out.close()
