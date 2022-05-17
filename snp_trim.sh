#!/bin/bash
#SBATCH --account=def-pawilson
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bcrosby@trentu.ca
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH --time=0-2:0:0
#SBATCH -o snp_trim.log
#SBATCH -e snp_trim.err


#
# This script requires the name of the raw fastq file directory as an argument,
# e.g, "sbatch snp_trim.sh MiSeq_12_microsatellite"
#


module load fastqc
module load java
module load trimmomatic


rm -fr $1_snp_trim/
mkdir $1_snp_trim/
mkdir $1_snp_trim/unpaired/


echo "########################################################"
echo "# Trimming data for run: $1 "
echo "########################################################"

echo "########################################################" >> snp_trim.err
echo "# Trimming data for run: $1 " >> snp_trim.err
echo "########################################################" >> snp_trim.err


ls /home/bcrosby/projects/def-pawilson/MiSeq_microsatellite/caribou/$1/fastq/*R1*.gz | \
	sed -r "s:_S[0-9]+_L001_R1_001.fastq.gz::" | \
        sed -r "s:/home/bcrosby/projects/def-pawilson/MiSeq_microsatellite/caribou/.*/fastq/::" \
	> $1_snp_trim/sample_list.txt


while IFS= read -r SAMPLE; do

	echo "# Trimming sample ${SAMPLE} #"
	echo "# Trimming sample ${SAMPLE} #" >> snp_trim.err

	java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE \
		-threads 4 \
		/home/bcrosby/projects/def-pawilson/MiSeq_microsatellite/caribou/$1/fastq/${SAMPLE}_S[0-9]*_L001_R1_001.fastq.gz \
		/home/bcrosby/projects/def-pawilson/MiSeq_microsatellite/caribou/$1/fastq/${SAMPLE}_S[0-9]*_L001_R2_001.fastq.gz \
		./$1_snp_trim/${SAMPLE}_R1_snp_trim.fastq.gz \
		./$1_snp_trim/unpaired/${SAMPLE}_R1_unpaired.fastq.gz \
		./$1_snp_trim/${SAMPLE}_R2_snp_trim.fastq.gz \
		./$1_snp_trim/unpaired/${SAMPLE}_R2_unpaired.fastq.gz \
		ILLUMINACLIP:adapter_primer.fa:2:30:10 \
		SLIDINGWINDOW:4:10 MINLEN:110


	echo "# Trimming for ${SAMPLE} complete #"
	echo "# Trimming for ${SAMPLE} complete #" >> snp_trim.err

done < $1_snp_trim/sample_list.txt


mv snp_trim.log $1_snp_trim/
mv snp_trim.err $1_snp_trim/

