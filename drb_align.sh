#!/bin/bash
#SBATCH --account=def-pawilson
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bcrosby@trentu.ca
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH --time=0-1:00:0
#SBATCH -o drb_align_%A.log
#SBATCH -e drb_align_%A.err


#
# This script has two requirements:
#
#	1. The name of the sequencing run as an input argument
#		e.g, "sbatch drb_align.sh BRMS010"
#
#       2. output from trim.sh for that same sequencing run
#               Make sure that the argument entered for Requirement 1 is the
#               same as the argument entered when running trim.sh
#


module load bowtie2
module load samtools


rm -fr $1_drb_align/
mkdir $1_drb_align/
mkdir $1_drb_align/dropped/


# read_group=$(grep '<Flowcell>' ~/projects/def-pawilson/caribou_MiSeq_project/$1/RunInfo.xml | sed -r "s:.*(\w{5})</Flowcell>:--rgID \1\.1:")


echo "########################################################"
echo "# Generating DRB alignments for run: $1 "
echo "########################################################"

echo "########################################################" >> drb_align_$SLURM_JOB_ID.err
echo "# Generating DRB alignments for run: $1 " >> drb_align_$SLURM_JOB_ID.err
echo "########################################################" >> drb_align_$SLURM_JOB_ID.err


ls /home/bcrosby/projects/def-pawilson/caribou_MiSeq_project/$1/fastq/*R1*.gz | \
        sed -r "s:_S[0-9]+_L001_R1_001.fastq.gz::" | \
        sed -r "s:/home/bcrosby/projects/def-pawilson/caribou_MiSeq_project/.*/fastq/::" \
        > $1_drb_align/sample_list.txt


while IFS= read -r SAMPLE; do


        echo "Aligning sample ${SAMPLE}"
        echo "Aligning sample ${SAMPLE}" >> drb_align_$SLURM_JOB_ID.err


#	bowtie2 "$(< $1_drb_align/bowtie_args.txt)"

	bowtie2 --end-to-end \
		-x references/drb_caribou_consensus \
		-1 ./$1_snp_trim/${SAMPLE}_R1_snp_trim.fastq.gz \
		-2 ./$1_snp_trim/${SAMPLE}_R2_snp_trim.fastq.gz \
		-S ./$1_drb_align/${SAMPLE}_drb.sam \
		--phred33 \
		--un-conc-gz ./$1_drb_align/dropped/${SAMPLE}_drb_dropped.sam \
		--rg-id \
		--rg SM: \
		--rg LB: \
		--threads 4


#                --rg-id ${read_group} \
#                --rg SM:${SAMPLE} \
#                --rg LB:${SAMPLE}-1 \



#	head ./$1_drb_align/${SAMPLE}_drb_temp.sam -n 4 > \
#		./$1_drb_align/${SAMPLE}_drb.sam


#       grep 'GATGGATCCTCTCTCTGCAGCACATTTCCT' ./$1_drb_align/${SAMPLE}_drb_temp.sam >> \
#               ./$1_drb_align/${SAMPLE}_drb.sam


        samtools view -bS -q 30 -@ 3 ./$1_drb_align/${SAMPLE}_drb.sam | \
                samtools sort -@ 4 - -o ./$1_drb_align/${SAMPLE}_drb.bam


	echo "Alignment for sample ${SAMPLE} complete"
	echo "Alignment for sample ${SAMPLE} complete" >> drb_align_$SLURM_JOB_ID.err


done < $1_drb_align/sample_list.txt


rm -r $1_drb_align/dropped/


mv drb_align_$SLURM_JOB_ID.log $1_drb_align/
mv drb_align_$SLURM_JOB_ID.err $1_drb_align/
