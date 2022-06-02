#!/bin/bash
#SBATCH --account=def-pawilson
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bcrosby@trentu.ca
#SBATCH -c 4
#SBATCH --mem=2G
#SBATCH --time=0-0:10:0
#SBATCH -o primer_search_%A.log
#SBATCH -e primer_search_%A.err


#
# This script requires the name of the raw fastq file directory as an argument,
# and the name of a forward primer as an argument
# e.g, "sbatch primer_search.sh BRMS010 BM4513"
#

rm -fr $1_primer_search_$2/
mkdir $1_primer_search_$2/


echo "########################################################"
echo "# Searching for $2 primer sequences in run: $1 "
echo "########################################################"

echo "########################################################" >> primer_search_$SLURM_JOB_ID.err
echo "# Searching for $2 primer sequences in run: $1 " >> primer_search_$SLURM_JOB_ID.err
echo "########################################################" >> primer_search_$SLURM_JOB_ID.err


ls /home/bcrosby/projects/def-pawilson/caribou_MiSeq_project/$1/fastq/*R1*.gz | \
	sed -r "s:_S[0-9]+_L001_R1_001.fastq.gz::" | \
        sed -r "s:/home/bcrosby/projects/def-pawilson/caribou_MiSeq_project/.*/fastq/::" \
	> $1_primer_search_$2/sample_list.txt


grep "$2" primers_f.txt | sed -r "s/$2\t([A-Z]+)/\1/" > temp_query_$2-$1.txt


while IFS= read -r SAMPLE; do


	echo "# Searching sample ${SAMPLE} #"
	echo "# Searching sample ${SAMPLE} #" >> primer_search_$SLURM_JOB_ID.err


	zcat  /home/bcrosby/projects/def-pawilson/caribou_MiSeq_project/$1/fastq/${SAMPLE}_S[0-9]*_L001_R1_001.fastq.gz | \
		grep -f temp_query_$2-$1.txt | \
		sort | \
		uniq -c | \
		sort -bgr > \
		./$1_primer_search_$2/${SAMPLE}_$2.txt


done < $1_primer_search_$2/sample_list.txt


rm temp_query_$2-$1.txt

mv primer_search_$SLURM_JOB_ID.log $1_primer_search_$2/
mv primer_search_$SLURM_JOB_ID.err $1_primer_search_$2/
