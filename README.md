# MiSeq_microsatellite
Repository for maintaining scripts and tools used in processing and analysis of caribou (Rangifer tarandus) high-throughput sequencing amplicon data


###  These scripts are used for processing and genotyping high-throughput amplicon sequences in caribou.


### trim.sh

Used for quality-trimming and quality assessment. **Trimmomatic** removes adapter sequences, removes low-quality sequences, and removes truncated sequences.  **Fastqc** generates quality reports for trimmed sequence data


### genome_align.sh

**Bowtie2** aligns trimmed sequence reads to the caribou reference genome. **Samtools** takes SAM output from Bowtie2 and converts it to a set of sorted, compressed BAM files.


### megasat.sh

megasat.sh runs **MEGASAT_Genotype.pl** (Zhan *et al,* 2017, https://github.com/beiko-lab/MEGASAT), taking the output from trim.sh as input and generating a file of microsatellite genotypes.


### sexID.sh

**Bowtie2** aligns trimmed sequence reads to the caribou ZFX/Y locus for sex identification. Bowtie2 SAM output is converted to BAM using **samtools*, then a custom script generates a sexID table for all samples run.
