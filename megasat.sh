#!/bin/bash
#SBATCH --account=def-pawilson
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bcrosby@trentu.ca
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH --time=0-2:0:0
#SBATCH -o megasat.log
#SBATCH -e megasat.err


#
# This script requires the name of the fastq file directory as an input argument,
# e.g, "sbatch megasat.sh MiSeq_12_microsatellite"
#


module load perl


rm -fr ./$1_merged/
mkdir ./$1_merged/

rm -fr ./$1_megasat/
mkdir ./$1_megasat/


echo "#############################################################"
echo "# Running MEGASAT_Genotype.pl on run: $1 "
echo "#############################################################"

echo "#############################################################" >> megasat.err
echo "# Running MEGASAT_Genotype.pl on run: $1 " >> megasat.err
echo "#############################################################" >> megasat.err


ls /home/bcrosby/projects/def-pawilson/MiSeq_microsatellite/caribou/$1/fastq/*R1*.gz | \
        sed -r "s:_S[0-9]+_L001_R1_001.fastq.gz::" | \
        sed -r "s:/home/bcrosby/projects/def-pawilson/MiSeq_microsatellite/caribou/.*/fastq/::" \
        > $1_megasat/sample_list.txt


while IFS= read -r SAMPLE; do

        ((i=i%4)); ((i++==0)) && wait

        echo "# Merging sample ${SAMPLE} #" &
        echo "# Merging sample ${SAMPLE} #" >> megasat.err &

	/home/bcrosby/projects/def-pawilson/software/usearch11.0.667_i86linux32 \
		-fastq_mergepairs ./$1_trim/${SAMPLE}_R1_trim.fastq -fastqout ./$1_merged/${SAMPLE}_merged.fastq &

done < $1_megasat/sample_list.txt


perl /home/bcrosby/projects/def-pawilson/software/MEGASAT-master/'MEGASAT_1.0 for Linux'/MEGASAT_Genotype.pl \
	./primerfile_1.2.txt \
	4 \
	50 \
	4 \
	./$1_merged/ \
	./$1_megasat/


mv megasat.log $1_megasat/
mv megasat.err $1_megasat/
