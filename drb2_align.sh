#!/bin/bash
#SBATCH --account=def-pawilson
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bcrosby@trentu.ca
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH --time=0-1:00:0
#SBATCH -o drb2_align.log
#SBATCH -e drb2_align.err


#
# This script has two requirements:
#
#	1. The name of the sequencing run as an input argument
#		e.g, "sbach drb2_align.sh MiSeq_12_microsatellite"
#
#       2. output from trim.sh for that same sequencing run
#               Make sure that the argument entered for Requirement 1 is the
#               same as the argument entered when running trim.sh
#


module load bowtie2
module load samtools


rm -fr drb2_align_$1/
mkdir drb2_align_$1/
mkdir drb2_align_$1/
mkdir drb2_align_$1/dropped/


echo "########################################################"
echo "# Generating DRB alignments for run: $1 "
echo "########################################################"

echo "########################################################" >> drb2_align.err
echo "# Generating DRB alignments for run: $1 " >> drb2_align.err
echo "########################################################" >> drb2_align.err


ls /home/bcrosby/projects/def-pawilson/MiSeq_microsatellite/caribou/$1/fastq/*R1*.gz | \
        sed -r "s:_S[0-9]+_L001_R1_001.fastq.gz::" | \
        sed -r "s:/home/bcrosby/projects/def-pawilson/MiSeq_microsatellite/caribou/.*/fastq/::" \
        > drb2_align_$1/sample_list.txt


while IFS= read -r SAMPLE; do


        echo "Aligning sample ${SAMPLE}"
        echo "Aligning sample ${SAMPLE}" >> drb2_align.err


        bowtie2 --end-to-end \
                -x references/drb_caribou_2 \
                -1 ./trim_$1/${SAMPLE}_R1_trim.fastq.gz \
                -2 ./trim_$1/${SAMPLE}_R2_trim.fastq.gz \
                -S ./drb2_align_$1/${SAMPLE}_drb.sam \
                --phred33 \
                --un-conc-gz ./drb2_align_$1/dropped/${SAMPLE}_drb_dropped.sam \
                --rg-id KBP3B.1 \
                --rg SM:${SAMPLE} \
                --rg LB:${SAMPLE}-1 \
		--threads 4


        samtools view -bS -q 30 -@ 3 ./drb2_align_$1/${SAMPLE}_drb.sam | \
                samtools sort -@ 4 - -o ./drb2_align_$1/${SAMPLE}_drb.bam


	echo "Alignment for sample ${SAMPLE} complete"
	echo "Alignment for sample ${SAMPLE} complete" >> drb2_align.err



done < drb2_align_$1/sample_list.txt


rm drb2_align_$1/read_lengths_temp.txt
rm -r drb2_align_$1/dropped/


mv drb2_align.log drb2_align_$1/
mv drb2_align.err drb2_align_$1/


rm -r drb2_align_$1/dropped/
