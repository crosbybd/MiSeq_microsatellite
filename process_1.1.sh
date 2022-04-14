#!/bin/bash
#SBATCH --account=def-pawilson
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bcrosby@trentu.ca
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH --time=0-5:0:0
#SBATCH -o process_1.1.log
#SBATCH -e process_1.1.err


module load fastqc
module load java
module load trimmomatic
module load bowtie2
module load samtools


rm -r trim_$1/
mkdir trim_$1/
mkdir trim_$1/unpaired/

rm -r sam_$1/
mkdir sam_$1/
mkdir sam_$1/dropped/

rm -r bam_$1/
mkdir bam_$1/

rm -r fastqc_$1/
mkdir fastqc_$1/


echo "########################################################"
echo "# Processing data for run: $1 "
echo "########################################################"

echo "########################################################" >> process_1.1.err
echo "# Processing data for run: $1 " >> process_1.1.err
echo "########################################################" >> process_1.1.err


ls /home/bcrosby/projects/def-pawilson/MiSeq_microsatellite/caribou/$1/fastq/*R1*.gz > fastq_list_R1.txt


sed -r "s:_S[0-9]+_L001_R1_001.fastq.gz::" fastq_list_R1.txt | \
        sed -r "s:/home/bcrosby/projects/def-pawilson/MiSeq_microsatellite/caribou/$1/fastq/::" | \
	> sample_list.txt


while IFS= read -r SAMPLE; do

	echo "# Processing sample ${SAMPLE} #"
	echo "# Processing sample ${SAMPLE} #" >> process_1.1.err

	java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE \
		-threads 4 \
		/home/bcrosby/projects/def-pawilson/MiSeq_microsatellite/caribou/$1/fastq/${SAMPLE}_S[0-9]+_L001_R1_001.fastq.gz \
		/home/bcrosby/projects/def-pawilson/MiSeq_microsatellite/caribou/$1/fastq/${SAMPLE}_S[0-9]+_L001_R2_001.fastq.gz \
		./trim_$1/${SAMPLE}_R1_trim.fastq.gz \
		./trim_$1/unpaired/${SAMPLE}_R1_unpaired.fastq.gz \
		./trim_$1/${SAMPLE}_R2_trim.fastq.gz \
		./trim_$1/unpaired/${SAMPLE}_R2_unpaired.fastq.gz \
		ILLUMINACLIP:adapter.fa:2:30:10 \
		SLIDINGWINDOW:4:10 MINLEN:110

	gunzip -k ./trim_$1/${SAMPLE}_R1_trim.fastq.gz &
	gunzip -k ./trim_$1/${SAMPLE}_R2_trim.fastq.gz

	fastqc ./trim_$1/${SAMPLE}_R1_trim.fastq.gz \
		-o ./fastqc_$1/

	bowtie2 --end-to-end \
		-x /home/bcrosby/projects/def-pawilson/reference_genome/Dovetail_hirise_May2021_final_assembly \
		-1 ./trim_$1/${SAMPLE}_R1_trim.fastq.gz \
		-2 ./trim_$1/${SAMPLE}_R2_trim.fastq.gz \
		-S ./sam_$1/${SAMPLE}.sam \
		--phred33 \
		--un-conc-gz ./sam_$1/dropped/${SAMPLE}_dropped.sam \
		--rg-id KBP3B.1 \
		--rg SM:${SAMPLE} \
		--rg LB:${SAMPLE}-1 \
		--threads 4

	samtools view -bS -q 30 -@ 3 ./sam_$1/${SAMPLE}.sam | \
		samtools sort - -@ 4 -o ./bam_$1/${SAMPLE}.bam

	echo "# Processing for ${SAMPLE} complete #"
	echo "# Processing for ${SAMPLE} complete #" >> process_1.1.err

done < sample_list.txt
