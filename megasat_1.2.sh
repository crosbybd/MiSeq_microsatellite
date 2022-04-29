#!/bin/bash
#SBATCH --account=def-pawilson
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bcrosby@trentu.ca
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH --time=0-1:0:0
#SBATCH -o megasat_1.2.log
#SBATCH -e megasat_1.2.err


#
# This script requires the name of the fastq file directory as an input argument,
# e.g, "sbatch megasat_1.2.sh MiSeq_12_microsatellite"
#


module load perl


rm -fr ./merged_$1/
mkdir ./merged_$1/

rm -fr ./megasat_$1/
mkdir ./megasat_$1/


echo "#############################################################"
echo "# Running MEGASAT_Genotype.pl on run: $1 "
echo "#############################################################"

echo "#############################################################" >> megasat_1.2.err
echo "# Running MEGASAT_Genotype.pl on run: $1 " >> megasat_1.2.err
echo "#############################################################" >> megasat_1.2.err


ls /home/bcrosby/projects/def-pawilson/MiSeq_microsatellite/caribou/$1/fastq/*R1*.gz | \
        sed -r "s:_S[0-9]+_L001_R1_001.fastq.gz::" | \
        sed -r "s:/home/bcrosby/projects/def-pawilson/MiSeq_microsatellite/caribou/.*/fastq/::" \
        > megasat_$1/sample_list.txt


while IFS= read -r SAMPLE; do

        ((i=i%4)); ((i++==0)) && wait

        echo "# Merging sample ${SAMPLE} #" &
        echo "# Merging sample ${SAMPLE} #" >> megasat_1.2.err &

	/home/bcrosby/projects/def-pawilson/software/usearch11.0.667_i86linux32 \
		-fastq_mergepairs ./trim/${SAMPLE}_R1_trim.fastq -fastqout ./merged_$1/${SAMPLE}_merged.fastq &

done < megasat_$1/sample_list.txt


perl /home/bcrosby/projects/def-pawilson/software/MEGASAT-master/'MEGASAT_1.0 for Linux'/MEGASAT_Genotype.pl \
	./primerfile_1.2.txt \
	4 \
	50 \
	4 \
	./merged_$1/ \
	./megasat_$1/


mv megasat_1.2.log megasat_$1/
mv megasat_1.2.err megasat_$1/
