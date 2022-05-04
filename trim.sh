#!/bin/bash
#SBATCH --account=def-pawilson
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bcrosby@trentu.ca
#SBATCH -c 4
#SBATCH --mem=8G
#SBATCH --time=0-2:0:0
#SBATCH -o trim.log
#SBATCH -e trim.err


#
# This script requires the name of the raw fastq file directory as an argument,
# e.g, "sbatch trim.sh MiSeq_12_microsatellite"
#


module load fastqc
module load java
module load trimmomatic


rm -fr trim_$1/
mkdir trim_$1/
mkdir trim_$1/unpaired/

rm -fr fastqc_$1/
mkdir fastqc_$1/


echo "########################################################"
echo "# Trimming data for run: $1 "
echo "########################################################"

echo "########################################################" >> trim.err
echo "# Trimming data for run: $1 " >> trim.err
echo "########################################################" >> trim.err


ls /home/bcrosby/projects/def-pawilson/MiSeq_microsatellite/caribou/$1/fastq/*R1*.gz | \
	sed -r "s:_S[0-9]+_L001_R1_001.fastq.gz::" | \
        sed -r "s:/home/bcrosby/projects/def-pawilson/MiSeq_microsatellite/caribou/.*/fastq/::" \
	> trim_$1/sample_list.txt


while IFS= read -r SAMPLE; do

	echo "# Trimming sample ${SAMPLE} #"
	echo "# Trimming sample ${SAMPLE} #" >> trim.err

	java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE \
		-threads 4 \
		/home/bcrosby/projects/def-pawilson/MiSeq_microsatellite/caribou/$1/fastq/${SAMPLE}_S[0-9]*_L001_R1_001.fastq.gz \
		/home/bcrosby/projects/def-pawilson/MiSeq_microsatellite/caribou/$1/fastq/${SAMPLE}_S[0-9]*_L001_R2_001.fastq.gz \
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


	echo "# Trimming for ${SAMPLE} complete #"
	echo "# Trimming for ${SAMPLE} complete #" >> trim.err

done < trim_$1/sample_list.txt

mv trim.log trim_$1/
mv trim.err trim_$1/
