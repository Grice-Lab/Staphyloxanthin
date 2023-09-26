# Amy Campbell
# Updated for Staphyloxanthin 2022 paper 
# In light of the fact that none of the 'followup' contigs were non-S. aureus
# This script adds the 'follow up contigs' back to the other cleaned contigs using cat

relativepathData="../Data/"
fullpathData=$(echo "$(cd "$(dirname "$relativepathData")"; pwd)/$(basename "$relativepathData")")

outputdir=$fullpathData"/WGSData/FinalIlluminaContigs/"

clean=$fullpathData"/WGSData/ShortAssemblyCleanContigs/"

followuppath=$fullpathData"/WGSData/ShortAssemblyFollowUpContigs/"


mkdir -p $outputdir
"
for f in $clean*_CleanContigs.fasta; do

	bass=$(basename $f) 

	cleanext="_CleanContigs.fasta" 

	followupext="_FollowUpContigs.fasta"

	blank=""

	noext=${bass/$cleanext/$blank}

	#noextfullpath=${f/$cleanext/$blank}	
	
	newext="_cleaned.fasta"

	newname=$noext$newext
	
	cat $f $followuppath$noext$followupext > $outputdir$newname

done
