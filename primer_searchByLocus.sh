#!/bin/bash
#SBATCH --account=def-pawilson
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bcrosby@trentu.ca
#SBATCH -c 4
#SBATCH --mem=2G
#SBATCH --time=0-0:10:0
#SBATCH -o primer_search_%A.log


#
# This script requires the name of the raw fastq file directory as an argument,
#
# This script has two requirements, both entered as arguments:
#	1. The name of the directory containing raw fastq files, and
#	2. The name of the locus the user wishes to search for
# e.g, "sbatch primer_search.sh BRMS010 BM4513"
#


# Evaluate whether the primer_search directory exists for this dataset (specified in Argument 1). If not, create it
[ ! -d "$1_primer_search/" ] && mkdir $1_primer_search/


# Create a subdirectory which will contain the output for the locus specified in Argument 2
rm -fr $1_primer_search/$2/
mkdir $1_primer_search/$2/


echo "########################################################"
echo "# Searching for $2 primer sequences in run: $1 "
echo "########################################################"


# List the files in the directory specified in Argument 1, and convert that list to a set of sample names
ls /home/bcrosby/projects/def-pawilson/caribou_MiSeq_project/$1/fastq/*R1*.gz | \
	sed -r "s:_S[0-9]+_L001_R1_001.fastq.gz::" | \
        sed -r "s:/home/bcrosby/projects/def-pawilson/caribou_MiSeq_project/.*/fastq/::" \
	> $1_primer_search/$2/sample_list.txt


# Searches the default primers file for a match to the locus name specified in Argument 2, then stores its forward primer sequence
grep "$2" primers_f.txt | sed -r "s/$2\t([A-Z]+)/\1/" > temp_query_$2-$1.txt


# Iterate over each line of the sample_list.txt file
# Each line contains a sample name, so this loop operates on each sample, one at a time
while IFS= read -r SAMPLE; do


	echo "# Searching sample ${SAMPLE} #"


	# This command prints a sample's forward raw fastq.gz file and outputs lines that match the primer sequence in the temp_query file
	# Then, unique sequences are counted, and sequence counts are sorted by number, in descending order
	zcat  /home/bcrosby/projects/def-pawilson/caribou_MiSeq_project/$1/fastq/${SAMPLE}_S[0-9]*_L001_R1_001.fastq.gz | \
		grep -f temp_query_$2-$1.txt | \
		sort | \
		uniq -c | \
		sort -bgr > \
		./$1_primer_search/$2/${SAMPLE}_$2.txt


done < $1_primer_search/$2/sample_list.txt


# Delete the file containing the query sequence
rm temp_query_$2-$1.txt


# Move the job log file into the locus-specific subdirectory for safe keeping
mv primer_search_$SLURM_JOB_ID.log $1_primer_search/$2/
