#!/bin/bash
#SBATCH --account=def-pawilson
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bcrosby@trentu.ca
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH --time=0-1:00:0
#SBATCH -o prnp_vcf_%A.log
#SBATCH -e prnp_vcf_%A.err


module load samtools
module load java
module load gatk/4.2.4.0


rm -fr $1_prnp_vcf/
mkdir $1_prnp_vcf/


echo "-R references/prnp_exon3.fa" > $1_prnp_vcf/CombineGVCFs_args.txt
echo "-O $1_prnp_vcf/combined.g.vcf.gz" >> $1_prnp_vcf/CombineGVCFs_args.txt


ls /home/bcrosby/projects/def-pawilson/caribou_MiSeq_project/$1/fastq/*R1*.gz | \
        sed -r "s:_S[0-9]+_L001_R1_001.fastq.gz::" | \
        sed -r "s:/home/bcrosby/projects/def-pawilson/caribou_MiSeq_project/.*/fastq/::" \
        > $1_prnp_vcf/sample_list.txt


while IFS= read -r SAMPLE; do

	echo "Creating GVCF for sample ${SAMPLE}"
	echo "Creating GVCF for sample ${SAMPLE}" >> prnp_vcf_$SLURM_JOB_ID.err


	samtools index -b -@ 3 $1_prnp_align/${SAMPLE}_prnp.bam

	gatk --java-options "-Xmx16g -XX:ParallelGCThreads=4" \
		HaplotypeCaller \
		-R ./references/prnp_exon3.fa \
		-I $1_prnp_align/${SAMPLE}_prnp.bam \
		-O $1_prnp_vcf/${SAMPLE}.g.vcf.gz \
		-ERC GVCF


	echo "--variant $1_prnp_vcf/${SAMPLE}.g.vcf.gz" >> $1_prnp_vcf/CombineGVCFs_args.txt


done < $1_prnp_vcf/sample_list.txt


gatk --java-options "-Xmx16g -XX:ParallelGCThreads=4" \
	CombineGVCFs \
	--arguments_file $1_prnp_vcf/CombineGVCFs_args.txt


gatk --java-options "-Xmx16g -XX:ParallelGCThreads=4" \
	GenotypeGVCFs \
	-R ./references/prnp_exon3.fa \
	-V $1_prnp_vcf/$1_combined.g.vcf.gz \
	-O $1_prnp_vcf/$1_genotyped.g.vcf.gz


mv prnp_vcf_$SLURM_JOB_ID.log $1_prnp_vcf/
mv prnp_vcf_$SLURM_JOB_ID.err $1_prnp_vcf/

