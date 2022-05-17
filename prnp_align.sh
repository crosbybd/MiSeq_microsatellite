#!/bin/bash
#SBATCH --account=def-pawilson
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bcrosby@trentu.ca
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH --time=0-1:00:0
#SBATCH -o prnp_align.log
#SBATCH -e prnp_align.err


#
# This script has two requirements:
#
#	1. The name of the sequencing run as an input argument
#		e.g, "sbach prnp_align.sh MiSeq_12_microsatellite"
#
#       2. output from trim.sh for that same sequencing run
#               Make sure that the argument entered for Requirement 1 is the
#               same as the argument entered when running trim.sh
#


module load bowtie2
module load samtools


rm -fr $1_prnp_align/
mkdir $1_prnp_align/
mkdir $1_prnp_align/
mkdir $1_prnp_align/dropped/


echo "########################################################"
echo "# Generating PRNP genotypes for run: $1 "
echo "########################################################"

echo "########################################################" >> prnp_align.err
echo "# Generating PRNP genotypes for run: $1 " >> prnp_align.err
echo "########################################################" >> prnp_align.err


ls /home/bcrosby/projects/def-pawilson/MiSeq_microsatellite/caribou/$1/fastq/*R1*.gz | \
        sed -r "s:_S[0-9]+_L001_R1_001.fastq.gz::" | \
        sed -r "s:/home/bcrosby/projects/def-pawilson/MiSeq_microsatellite/caribou/.*/fastq/::" \
        > $1_prnp_align/sample_list.txt


while IFS= read -r SAMPLE; do


        echo "Aligning sample ${SAMPLE}"
        echo "Aligning sample ${SAMPLE}" >> prnp_align.err


        bowtie2 --end-to-end \
                -x references/prnp_exon3 \
                -1 ./$1_trim/${SAMPLE}_R1_trim.fastq.gz \
                -2 ./$1_trim/${SAMPLE}_R2_trim.fastq.gz \
                -S ./$1_prnp_align/${SAMPLE}_prnp.sam \
                --phred33 \
                --un-conc-gz ./$1_prnp_align/dropped/${SAMPLE}_prnp_dropped.sam \
                --rg-id KBP3B.1 \
                --rg SM:${SAMPLE} \
                --rg LB:${SAMPLE}-1 \
		--threads 4


        samtools view -bS -q 30 -@ 3 ./$1_prnp_align/${SAMPLE}_prnp.sam | \
                samtools sort -@ 4 - -o ./$1_prnp_align/${SAMPLE}_prnp.bam


	echo "Alignment for sample ${SAMPLE} complete"
	echo "Alignment for sample ${SAMPLE} complete" >> prnp_align.err



done < $1_prnp_align/sample_list.txt


rm $1_prnp_align/read_lengths_temp.txt
rm -r $1_prnp_align/dropped/


mv prnp_align.log $1_prnp_align/
mv prnp_align.err $1_prnp_align/


rm -r $1_prnp_align/dropped/
