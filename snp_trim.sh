#!/bin/bash
#SBATCH --account=def-pawilson
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bcrosby@trentu.ca
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH --time=0-2:0:0
#SBATCH -o snp_trim_%A.log


#
# This script requires the name of the raw fastq file directory as an argument,
# e.g, "sbatch snp_trim.sh BRMS010"
#


# Load the necessary functions into the environment

# Java is necessary for Trimmomatic
module load java

# Trimmomatic to perform quality-trimming of sequence reads
module load trimmomatic


# Create directory to contain trimmed sequences
rm -fr $1_snp_trim/
mkdir $1_snp_trim/
mkdir $1_snp_trim/unpaired/


echo "########################################################"
echo "# Trimming data for run: $1 "
echo "########################################################"


# List the files in the specified directory, then convert those filenames to sample names
ls /home/bcrosby/projects/def-pawilson/caribou_MiSeq_project/$1/fastq/*R1*.gz | \
	sed -r "s:_S[0-9]+_L001_R1_001.fastq.gz::" | \
        sed -r "s:/home/bcrosby/projects/def-pawilson/caribou_MiSeq_project/.*/fastq/::" \
	> $1_snp_trim/sample_list.txt


# Iterate over each sample name found in sample_list.txt, performing trimming and generating quality reports one-by-one for each sample
while IFS= read -r SAMPLE; do


	echo "# Trimming sample ${SAMPLE} #"


	# Trimmomatic to quality-trim sequence reads
	java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE \	# Operate on paired-end sequence data
		-threads 4 \	# Use 4 cores to speed up process
		/home/bcrosby/projects/def-pawilson/caribou_MiSeq_project/$1/fastq/${SAMPLE}_S[0-9]*_L001_R1_001.fastq.gz \	# Call the file containing raw forward reads
		/home/bcrosby/projects/def-pawilson/caribou_MiSeq_project/$1/fastq/${SAMPLE}_S[0-9]*_L001_R2_001.fastq.gz \	# Call the file containing raw reverse reads
		./$1_snp_trim/${SAMPLE}_R1_snp_trim.fastq.gz \	# Define the filename containing trimmed forward reads
		./$1_snp_trim/unpaired/${SAMPLE}_R1_unpaired.fastq.gz \	# Define the filename containing dropped forward reads
		./$1_snp_trim/${SAMPLE}_R2_snp_trim.fastq.gz \	# Define the filename containing trimmed reverse reads
		./$1_snp_trim/unpaired/${SAMPLE}_R2_unpaired.fastq.gz \	# Define the filename containing dropped reverse reads
		ILLUMINACLIP:adapter_primer.fa:2:30:10 \	# Remove sequences found in adapter.fa, set quality thresholds
		SLIDINGWINDOW:4:10 MINLEN:110	# Use sliding window to measure quality of bases, remove any sequence reads less than 110 bp


	echo "# Trimming for ${SAMPLE} complete #"


done < $1_snp_trim/sample_list.txt


# Move the job log file ito the snp_trim directory for safe keeping
mv snp_trim_$SLURM_JOB_ID.log $1_snp_trim/

