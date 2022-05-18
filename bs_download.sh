#!/bin/bash
#SBATCH --account=def-pawilson
#SBATCH --mail-type=ALL
#SBATCH --mail-user=bcrosby@trentu.ca
#SBATCH -c 2
#SBATCH --mem=4G
#SBATCH --time=0-1:0:0
#SBATCH -o bs_download.log
#SBATCH -e bs_download.err


#
# This script requires the name of the raw fastq file directory as an argument,
# e.g, "sbatch bs_download.sh BRMS010"
#



rm -fr /home/bcrosby/projects/def-pawilson/caribou_MiSeq_project/$1/
mkdir /home/bcrosby/projects/def-pawilson/caribou_MiSeq_project/$1/
mkdir /home/bcrosby/projects/def-pawilson/caribou_MiSeq_project/$1/fastq/


rm -fr $1_bs_download/
mkdir $1_bs_download/



echo "########################################################"
echo "# Downloading data from BaseSpace for run: $1 "
echo "########################################################"

echo "########################################################" >> bs_download.err
echo "# Downloading data from BaseSpace for run: $1 "  >> bs_download.err
echo "########################################################" >> bs_download.err



/home/bcrosby/projects/def-pawilson/software/bs download project \
	-n $1 \
	-q \
	-o /home/bcrosby/projects/def-pawilson/caribou_MiSeq_project/$1/fastq/


mv /home/bcrosby/projects/def-pawilson/caribou_MiSeq_project/$1/fastq/*/*.fastq.gz \
	/home/bcrosby/projects/def-pawilson/caribou_MiSeq_project/$1/fastq/


rm -fr /home/bcrosby/projects/def-pawilson/caribou_MiSeq_project/$1/fastq/*ds*/


/home/bcrosby/projects/def-pawilson/software/bs download run \
	-n $1 \
	-q \
	-o /home/bcrosby/projects/def-pawilson/caribou_MiSeq_project/$1/ \
	--extension="csv"


/home/bcrosby/projects/def-pawilson/software/bs download run \
	-n $1 \
	-q \
	-o /home/bcrosby/projects/def-pawilson/caribou_MiSeq_project/$1/ \
	--exclude="*" \
	--include="RunInfo.xml"


echo "########################################################"
echo "# BaseSpace download complete for run: $1 "
echo "########################################################"

echo "########################################################" >> bs_download.err
echo "# BaseSpace download complete for run: $1 " >> bs_download.err
echo "########################################################" >> bs_download.err


mv bs_download.log $1_bs_download/
mv bs_download.err $1_bs_download/

