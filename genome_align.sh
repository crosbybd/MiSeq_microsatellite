#!/bin/bash
#SBATCH --account=def-pawilson
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bcrosby@trentu.ca
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH --time=0-2:0:0
#SBATCH -o genome_align.log
#SBATCH -e genome_align.err


#
# This script has two requirements:
#
# 	1. the name of the sequencing run as an input argument
#		 e.g, "sbatch genome_align.sh MiSeq_12_microsatellite"
#
#	2. output from trim_1.1.sh for that same sequencing run
#		Make sure that the argument entered for Requirement 1 is the
#		same as the argument entered when running trim_1.1.sh
#


module load bowtie2
module load samtools


rm -fr $1_genome_align/
mkdir $1_genome_align/
mkdir $1_genome_align/dropped/


echo "########################################################"
echo "# Aligning data for run: $1 "
echo "########################################################"

echo "########################################################" >> genome_align.err
echo "# Aligning data for run: $1 " >> genome_align.err
echo "########################################################" >> genome_align.err


rgID=$(grep '<Flowcell>' ~/projects/def-pawilson/caribou_MiSeq_project/$1/RunInfo.xml | sed -r "s:.*(\w{5})</Flowcell>:\1\.1:" | dos2unix)


ls /home/bcrosby/projects/def-pawilson/caribou_MiSeq_project/$1/fastq/*R1*.gz | \
        sed -r "s:_S[0-9]+_L001_R1_001.fastq.gz::" | \
        sed -r "s:/home/bcrosby/projects/def-pawilson/caribou_MiSeq_project/.*/fastq/::" \
        > $1_genome_align/sample_list.txt


while IFS= read -r SAMPLE; do

	echo "# Aligning sample ${SAMPLE} to reference genome #"
	echo "# Aligning sample ${SAMPLE} to reference genome #" >> genome_align.err

	bowtie2 --end-to-end \
		-x /home/bcrosby/projects/def-pawilson/reference_genome/Dovetail_hirise_May2021_final_assembly \
		-1 ./$1_trim/${SAMPLE}_R1_trim.fastq.gz \
		-2 ./$1_trim/${SAMPLE}_R2_trim.fastq.gz \
		-S ./$1_genome_align/${SAMPLE}.sam \
		--phred33 \
		--un-conc-gz ./$1_genome_align/dropped/${SAMPLE}_dropped.sam \
		--rg-id $rgID \
		--rg SM:${SAMPLE} \
		--rg LB:${SAMPLE}-1 \
		--threads 4

	samtools view -bS -q 30 -@ 3 ./$1_genome_align/${SAMPLE}.sam | \
		samtools sort - -@ 4 -o ./$1_genome_align/${SAMPLE}.bam

	echo "# Genome alignment for ${SAMPLE} complete #"
	echo "# Genome alignment for ${SAMPLE} complete #" >> genome_align.err

done < $1_genome_align/sample_list.txt


mv genome_align.log $1_genome_align/
mv genome_align.err $1_genome_align/
