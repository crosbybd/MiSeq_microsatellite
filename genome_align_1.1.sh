#!/bin/bash
#SBATCH --account=def-pawilson
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bcrosby@trentu.ca
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH --time=0-2:0:0
#SBATCH -o genome_align_1.1.log
#SBATCH -e genome_align_1.1.err


#
# This script has two requirements:
#
# 	1. the name of the sequencing run as an input argument
#		 e.g, "sbatch genome_align_1.1.sh MiSeq_12_microsatellite"
#
#	2. output from trim_1.1.sh for that same sequencing run
#		Make sure that the argument entered for Requirement 1 is the
#		same as the argument entered when running trim_1.1.sh
#


module load bowtie2
module load samtools


rm -fr alignments_$1/
mkdir alignments_$1/
mkdir alignments_$1/dropped/


echo "########################################################"
echo "# Aligning data for run: $1 "
echo "########################################################"

echo "########################################################" >> genome_align_1.1.err
echo "# Aligning data for run: $1 " >> genome_align_1.1.err
echo "########################################################" >> genome_align_1.1.err


ls /home/bcrosby/projects/def-pawilson/MiSeq_microsatellite/caribou/$1/fastq/*R1*.gz | \
        sed -r "s:_S[0-9]+_L001_R1_001.fastq.gz::" | \
        sed -r "s:/home/bcrosby/projects/def-pawilson/MiSeq_microsatellite/caribou/.*/fastq/::" \
        > alignments_$1/sample_list.txt


while IFS= read -r SAMPLE; do

	echo "# Aligning sample ${SAMPLE} to reference genome #"
	echo "# Aligning sample ${SAMPLE} to reference genome #" >> genome_align_1.1.err

	bowtie2 --end-to-end \
		-x /home/bcrosby/projects/def-pawilson/reference_genome/Dovetail_hirise_May2021_final_assembly \
		-1 ./trim_$1/${SAMPLE}_R1_trim.fastq.gz \
		-2 ./trim_$1/${SAMPLE}_R2_trim.fastq.gz \
		-S ./alignments_$1/${SAMPLE}.sam \
		--phred33 \
		--un-conc-gz ./alignments_$1/dropped/${SAMPLE}_dropped.sam \
		--rg-id KBP3B.1 \
		--rg SM:${SAMPLE} \
		--rg LB:${SAMPLE}-1 \
		--threads 4

	samtools view -bS -q 30 -@ 3 ./alignments_$1/${SAMPLE}.sam | \
		samtools sort - -@ 4 -o ./alignments_$1/${SAMPLE}.bam

	echo "# Genome alignment for ${SAMPLE} complete #"
	echo "# Genome alignment for ${SAMPLE} complete #" >> genome_align_1.1.err

done < alignments_$1/sample_list.txt


mv genome_align_1.1.log alignments_$1/
mv genome_align_1.1.err alignments_$1/
