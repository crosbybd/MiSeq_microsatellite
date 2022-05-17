#!/bin/bash
#SBATCH --account=def-pawilson
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bcrosby@trentu.ca
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH --time=0-1:00:0
#SBATCH -o drb_vcf.log
#SBATCH -e drb_vcf.err


module load samtools
module load java
module load gatk/4.2.4.0


rm -fr $1_drb_vcf/
mkdir $1_drb_vcf/


echo "-R references/drb_caribou.fasta" > $1_drb_vcf/CombineGVCFs_args.txt
echo "-O $1_drb_vcf/combined.g.vcf.gz" >> $1_drb_vcf/CombineGVCFs_args.txt


ls /home/bcrosby/projects/def-pawilson/MiSeq_microsatellite/caribou/$1/fastq/*R1*.gz | \
        sed -r "s:_S[0-9]+_L001_R1_001.fastq.gz::" | \
        sed -r "s:/home/bcrosby/projects/def-pawilson/MiSeq_microsatellite/caribou/.*/fastq/::" \
        > $1_drb_vcf/sample_list.txt


while IFS= read -r SAMPLE; do

	echo "Creating GVCF for sample ${SAMPLE}"
	echo "Creating GVCF for sample ${SAMPLE}" >> drb_vcf.err


	samtools index -b -@ 3 drb_align_$1/${SAMPLE}_drb.bam

	gatk --java-options "-Xmx16g -XX:ParallelGCThreads=4" \
		HaplotypeCaller \
		-R ./references/drb_caribou.fasta \
		-I drb_align_$1/${SAMPLE}_drb.bam \
		-O $1_drb_vcf/${SAMPLE}.g.vcf.gz \
		-ERC GVCF


	echo "--variant $1_drb_vcf/${SAMPLE}.g.vcf.gz" >> $1_drb_vcf/CombineGVCFs_args.txt


done < $1_drb_vcf/sample_list.txt


gatk --java-options "-Xmx16g -XX:ParallelGCThreads=4" \
	CombineGVCFs \
	--arguments_file $1_drb_vcf/CombineGVCFs_args.txt


gatk --java-options "-Xmx16g -XX:ParallelGCThreads=4" \
	GenotypeGVCFs \
	-R ./references/drb_caribou.fasta \
	-V $1_drb_vcf/combined.g.vcf.gz \
	-O $1_drb_vcf/genotyped.g.vcf.gz


mv drb_vcf.log $1_drb_vcf/
mv drb_vcf.err $1_drb_vcf/

