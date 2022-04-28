#!/bin/bash
#SBATCH --account=def-pawilson
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bcrosby@trentu.ca
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH --time=0-1:00:0
#SBATCH -o sexID_1.1.log
#SBATCH -e sexID_1.1.err


#
# This script has two requirements:
#
#	1. The name of the sequencing run as an input argument
#		e.g, "sbach sexID_1.1.sh MiSeq_12_microsatellite"
#
#       2. output from trim_1.1.sh for that same sequencing run
#               Make sure that the argument entered for Requirement 1 is the
#               same as the argument entered when running trim_1.1.sh
#


module load samtools


rm -r sexID_$1/
mkdir sexID_$1/
mkdir sexID_$1/alignments/
mkdir sexID_$1/alignments/dropped/


echo "########################################################"
echo "# Generating sexID table for run: $1 "
echo "########################################################"

echo "########################################################" >> sexID_1.1.err
echo "# Generating sexID table for run: $1 " >> sexID_1.1.err
echo "########################################################" >> sexID_1.1.err


ls /home/bcrosby/projects/def-pawilson/MiSeq_microsatellite/caribou/$1/fastq/*R1*.gz > fastq_list_R1.txt


sed -r "s:_S[0-9]+_L001_R1_001.fastq.gz::" fastq_list_R1.txt | \
        sed -r "s:/home/bcrosby/projects/def-pawilson/MiSeq_microsatellite/caribou/$1/fastq/::" | \
        > sample_list.txt


echo "sample,X_depth,Y_depth,sex" > sexID_$1/sexID_results_$1.csv


while IFS= read -r SAMPLE; do


        echo "Aligning sample ${SAMPLE}"
        echo "Aligning sample ${SAMPLE}" >> sexID_1.1.err


        bowtie2 --local \
                -x ./references/zfy_caribou \
                -1 ./trim_$1/${SAMPLE}_R1_trim.fastq.gz \
                -2 ./trim_$1/${SAMPLE}_R2_trim.fastq.gz \
                -S ./sexID_$1/alignments/${SAMPLE}_sexID.sam \
                --phred33 \
                --un-conc-gz ./sexID_$1/alignments/dropped/${SAMPLE}_sexID_dropped.sam \
                --rg-id KBP3B.1 \
                --threads 4


        samtools view -bS -q 30 -@ 3 ./sexID_$1/alignments/${SAMPLE}_sexID.sam | \
                samtools sort -@ 4 - -o ./sexID_$1/alignments/${SAMPLE}_sexID.bam


	echo "Alignment for sample ${SAMPLE} complete"
	echo "Alignment for sample ${SAMPLE} complete" >> sexID_1.1.err


        echo "Genotyping sample ${SAMPLE}"
        echo "Genotyping sample ${SAMPLE}" >> sexID_1.1.err


        samtools stats --reference references/zfy_caribou.fasta \
                -@ 3 \
                sexID_$1/alignments/${SAMPLE}_sexID.bam | \
                grep -E '^RL' | \
                cut -f 2,3 > \
                sexID_$1/read_lengths_temp.txt

                X_depth=$(grep '204' sexID_$1/read_lengths_temp.txt | \
                        cut -f 2)

                Y_depth=$(grep '185' sexID_$1/read_lengths_temp.txt | \
                        cut -f 2)


                if [[ $(( $X_depth / 10 )) -gt $Y_depth || -z "$Y_depth" ]]
                then
                        echo "${SAMPLE},$X_depth,$Y_depth,F" >> sexID_$1/sexID_results_$1.csv
                else
                        echo "${SAMPLE},$X_depth,$Y_depth,M" >> sexID_$1/sexID_results_$1.csv
                fi




done < sample_list.txt


rm sexID_$1/read_lengths_temp.txt
rm sample_list.txt
rm fastq_list_R1.txt
rm -r sexID_$1/alignments/dropped/
