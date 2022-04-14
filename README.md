# MiSeq_microsatellite
For maintaining scripts and tools used in processing and analysis of caribou (Rangifer tarandus) microsatellite amplicon data


###  The scripts herein are used for processing and genotyping caribou microsatellite amplicon data.


### process_1.1.sh

The primary function of this script is to quality-trim sequence data and align it to the reference genome. Here are the tools used in process_1.1.sh and their functions:


<ol>
	<li>**trimmomatic**		Remove adapter sequences, remove low-quality sequences, and remove truncated sequences</li>
	<li>**fastqc**		Generate quality reports for trimmed sequence data</li>
	<li>**bowtie2**		Align trimmed sequences to the caribou reference genome</li>
	<li>**samtools view**	Convert alignment output from bowtie2 from uncompressed SAM to compressed BAM and remove low-quality alignments.</li>
	<li>**samtools sort**	Sort sequence reads in the BAM file; typically required for downstream analysis.</li>
</ol>

### megasat_1.2.sh

megasat_1.2.sh runs MEGASAT_Genotype.pl (Zhan *et al,* 2017, https://github.com/beiko-lab/MEGASAT), taking the output from process_1.1.sh as input and generating a file of microsatellite genotypes.
