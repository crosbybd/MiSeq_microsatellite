#!/bin/bash
#SBATCH --account=def-pawilson
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bcrosby@trentu.ca
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH --time=0-1:00:0
#SBATCH -o sexID_%A.log
#SBATCH -e sexID_%A.err


#
# This script has two requirements:
#
#	1. The name of the sequencing run as an input argument
#		e.g, "sbach sexID.sh MiSeq_12_microsatellite"
#
#       2. output from trim_1.1.sh for that same sequencing run
#               Make sure that the argument entered for Requirement 1 is the
#               same as the argument entered when running trim_1.1.sh
#


module load bowtie2
module load samtools


rm -fr $1_sexID/
mkdir $1_sexID/
mkdir $1_sexID/alignments/
mkdir $1_sexID/alignments/dropped/


echo "########################################################"
echo "# Generating sexID table for run: $1 "
echo "########################################################"

echo "########################################################" >> sexID_$SLURM_JOB_ID.err
echo "# Generating sexID table for run: $1 " >> sexID_$SLURM_JOB_ID.err
echo "########################################################" >> sexID_$SLURM_JOB_ID.err


ls /home/bcrosby/projects/def-pawilson/caribou_MiSeq_project/$1/fastq/*R1*.gz | \
        sed -r "s:_S[0-9]+_L001_R1_001.fastq.gz::" | \
        sed -r "s:/home/bcrosby/projects/def-pawilson/caribou_MiSeq_project/.*/fastq/::" \
        > $1_sexID/sample_list.txt


echo "sample,X_depth,Y_depth,sex" > $1_sexID/sexID_results_$1.csv


while IFS= read -r SAMPLE; do


        echo "Aligning sample ${SAMPLE}"
        echo "Aligning sample ${SAMPLE}" >> sexID_$SLURM_JOB_ID.err


        bowtie2 --local \
                -x ./references/zfy_caribou \
                -1 ./$1_trim/${SAMPLE}_R1_trim.fastq.gz \
                -2 ./$1_trim/${SAMPLE}_R2_trim.fastq.gz \
                -S ./$1_sexID/alignments/${SAMPLE}_sexID.sam \
                --phred33 \
                --un-conc-gz ./$1_sexID/alignments/dropped/${SAMPLE}_sexID_dropped.sam \
                --rg-id KBP3B.1 \
                --threads 4


        samtools view -bS -q 30 -@ 3 ./$1_sexID/alignments/${SAMPLE}_sexID.sam | \
                samtools sort -@ 4 - -o ./$1_sexID/alignments/${SAMPLE}_sexID.bam


	echo "Alignment for sample ${SAMPLE} complete"
	echo "Alignment for sample ${SAMPLE} complete" >> sexID_$SLURM_JOB_ID.err


        echo "Genotyping sample ${SAMPLE}"
        echo "Genotyping sample ${SAMPLE}" >> sexID_$SLURM_JOB_ID.err


        samtools stats --reference references/zfy_caribou.fasta \
                -@ 3 \
                $1_sexID/alignments/${SAMPLE}_sexID.bam | \
                grep -E '^RL' | \
                cut -f 2,3 > \
                $1_sexID/read_lengths_temp.txt


		total_depth=$(cut -f 2 $1_sexID/read_lengths_temp.txt | awk '{sum += $1 } END { print sum }')

                X_depth=$(grep '204' $1_sexID/read_lengths_temp.txt | \
                        cut -f 2)

                Y_depth=$(grep '185' $1_sexID/read_lengths_temp.txt | \
                        cut -f 2)


                if [[ -z "$X_depth" || $total_depth -lt 20  ]]
                then

                        echo "${SAMPLE},$X_depth,$Y_depth,-99" >> $1_sexID/sexID_results_$1.csv

                elif [[ $(( $X_depth / 10 )) -gt $Y_depth || -z "$Y_depth" ]]
                then

                        echo "${SAMPLE},$X_depth,$Y_depth,F" >> $1_sexID/sexID_results_$1.csv

                else

                        echo "${SAMPLE},$X_depth,$Y_depth,M" >> $1_sexID/sexID_results_$1.csv

                fi

done < $1_sexID/sample_list.txt


rm $1_sexID/read_lengths_temp.txt
rm -r $1_sexID/alignments/dropped/


mv sexID_$SLURM_JOB_ID.log $1_sexID/
mv sexID_$SLURM_JOB_ID.err $1_sexID/
