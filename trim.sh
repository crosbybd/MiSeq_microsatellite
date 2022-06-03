#!/bin/bash
#SBATCH --account=def-pawilson
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bcrosby@trentu.ca
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH --time=0-2:0:0
#SBATCH -o trim_%A.log


#
# This script requires the name of the raw fastq file directory as an argument,
# e.g, "sbatch trim.sh BRMS010"
#


# Load the necessary functions into the environment

# Fastqc for generating sequence quality reports
module load fastqc

# Java is needed by Trimmomatic
module load java

# Trimmomatic for quality-filtering and trimming artificial reads from data
module load trimmomatic


# Create a directory to contain trimmed sequences for this dataset
rm -fr $1_trim/
mkdir $1_trim/
mkdir $1_trim/unpaired/


# Create a directory to contain the quality reports for this dataset
rm -fr $1_fastqc/
mkdir $1_fastqc/


echo "########################################################"
echo "# Trimming data for run: $1 "
echo "########################################################"


# List the files in the specified directory, then convert those filenames into sample names
ls /home/bcrosby/projects/def-pawilson/caribou_MiSeq_project/$1/fastq/*R1*.gz | \
	sed -r "s:_S[0-9]+_L001_R1_001.fastq.gz::" | \
        sed -r "s:/home/bcrosby/projects/def-pawilson/caribou_MiSeq_project/.*/fastq/::" \
	> $1_trim/sample_list.txt


# Iterate over each sample name in sample_list.txt, conducting quality-trimming and generating quality reports one by one for each sample
while IFS= read -r SAMPLE; do


	echo "# Trimming sample ${SAMPLE} #"


	# Call Trimmomatic
	java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE \	# Operate on paired-end sequences
		-threads 4 \	# Use 4 cores to speed up the process
		/home/bcrosby/projects/def-pawilson/caribou_MiSeq_project/$1/fastq/${SAMPLE}_S[0-9]*_L001_R1_001.fastq.gz \	# Call the input forward reads
		/home/bcrosby/projects/def-pawilson/caribou_MiSeq_project/$1/fastq/${SAMPLE}_S[0-9]*_L001_R2_001.fastq.gz \	# Call the input reverse reads
		./$1_trim/${SAMPLE}_R1_trim.fastq.gz \	# Define the name of the output file containing trimmed forward reads
		./$1_trim/unpaired/${SAMPLE}_R1_unpaired.fastq.gz \	# Define the name of the file containing forward reads that are dropped due to having low quality
		./$1_trim/${SAMPLE}_R2_trim.fastq.gz \	# Define the name of the output file containing trimmed reverse reads
		./$1_trim/unpaired/${SAMPLE}_R2_unpaired.fastq.gz \	# Define the name of the file containing reverse reads that are dropped due to having low quality
		ILLUMINACLIP:adapter.fa:2:30:10 \ # Remove sequences found in the adapter.fa file, and set thresholds for quality
		SLIDINGWINDOW:4:10 MINLEN:110	# Use a sliding window to assess quality and drop all reads shorter than 110 bp


	# Decompress trimmed sequence data
	gunzip -k ./$1_trim/${SAMPLE}_R1_trim.fastq.gz &
	gunzip -k ./$1_trim/${SAMPLE}_R2_trim.fastq.gz


	# Generate fastqc quality report for the forward read (read depth should be the same for forward and reverse reads)
	fastqc ./$1_trim/${SAMPLE}_R1_trim.fastq.gz \
		-o ./$1_fastqc/


	echo "# Trimming for ${SAMPLE} complete #"


done < $1_trim/sample_list.txt


# Activate virtual Python environment to access Multiqc
source ENV_multiqc/bin/activate


# Call Multiqc to compile all fastqc reports into a single document, put this document into the fastqc folder
multiqc ./$1_fastqc/ \
	-o ./$1_fastqc/$1_multiqc_report.html \
	-f


# Move job log file into trim folder for safe keeping
mv trim_$SLURM_JOB_ID.log $1_trim/
