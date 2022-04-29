# MiSeq_microsatellite
Repository for maintaining scripts and tools used in processing and analysis of caribou (Rangifer tarandus) high-throughput sequencing amplicon data


###  The scripts herein are used for processing and genotyping caribou microsatellite amplicon data.


### trim_1.1.sh

Used for quality-trimming and quality assessment. **Trimmomatic** removes adapter sequences, removes low-quality sequences, and removes truncated sequences.  **Fastqc** generates quality reports for trimmed sequence data


### genome_align_1.1.sh

**Bowtie2** aligns trimmed sequence reads to the caribou reference genome. **Samtools** takes SAM output from Bowtie2 and converts it to a set of sorted, compressed BAM files.


### megasat_1.2.sh

megasat_1.2.sh runs MEGASAT_Genotype.pl (Zhan *et al,* 2017, https://github.com/beiko-lab/MEGASAT), taking the output from process_1.1.sh as input and generating a file of microsatellite genotypes.


### sexID_1.1.sh

**Bowtie2** aligns trimmed sequence reads to the caribou ZFX/Y locus for sex identification. Bowtie2 SAM output is converted to BAM using **samtools*, then a custom script generates a sexID table for all samples run.
