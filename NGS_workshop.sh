##Alan Pittman December 2019

##################################         PART1      ###################################################       
###############Preparing our mitochondrial reference sequence fasta files################################
#########################################################################################################

## Download the human mitochondrial refernce sequence:
mkdir reference
wget --timestamping 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes/chrM.fa.gz' -O reference/chrM.fa.gz ###HAVE A BACKUP PLAN

## unzip the compressed chrM fasta txt file
gunzip reference/chrM.fa.gz

##view the top 10 lines of the file with the linux 'head' commands
head reference/chrM.fa

##Prepare th FASTA file for use as a referenence for BWA###

"""For the aligment algorithms of BWA, we first need to construct the FM-index for the reference genome 
(the index command)."""

software/bwa index reference/chrM.fa

##Prepare the FASTA file for use as a referenence for the genome analysis toolkit (GATK)###

"""Creat two needed files to access and safety check access to the reference files: a .dict dictionary 
of the contig names and sizes and a .fai fasta index file to allow efficient random access to the 
reference bases. You have to generate these files in order to be able to use a Fasta file as reference."""

##Creating the fasta index file

"""We use the faidx command in samtools to prepare the fasta index file. 
This file describes byte offsets in the fasta file for each contig, allowing us to compute exactly where 
a particular reference base at contig:pos is in the fasta file."""

software/samtools faidx reference/chrM.fa

"""This produces a text file with one record per line for each of the fasta contigs. Each record is of the: 
contig, size, location, basesPerLine, bytesPerLine."""

##Creating the sequence dictionary file

"""We use CreateSequenceDictionary.jar from Picard tools to create a .dict file from a fasta file."""

java -jar software/gatk-package-4.0.4.0-local.jar CreateSequenceDictionary -R reference/chrM.fa -O reference/chrM.dict

"""This produces a SAM-style header file; simply describing the contents of our fasta file."""

##################################         PART2      ###################################################       
################Quality control of our mitochondrial NGS sequence data prior to alignments###############
#########################################################################################################

##Trimming by quality and adapter removal with fastp

"""The pile of reads coming out of the sequencer is stored in FASTQ files. As a typical paired-end
sequencing experiment, each sequenced sample will produce two FASTQ files; One for the first pair of reads
and one for the second pair of reads."""

"""FASTQ files are usually quite big files containing information for all the reads produced. To get an insight
on how a FASTQ file look like we can simply open its first 12 lines:"""
head -n 12 sample1_r1.fastq 

"""With the program fastp we are assesing the quality of our sequence fastq data, which we can visualise in 
the output HTML report file"""

"""In addition to assesing the quality of our sequence data, fastp will trim off low quality trailing bases 
from our sequence reads and detect and remove and remaining unwanted adapter sequencenes"""

"""new cleaned fastq files are produced ready for the next steps"""
software/fastp -i sample1_r1.fastq -I sample1_r2.fastq -o sample1_out.R1.fq.gz -O sample1_out.R2.fq.gz --html \
 sample1_results.html --json sample1_results.json --report_title sample1_results


##Alternative way of visualising the quality control of the reads
"""An alternative way of visualising the quality of our sequence fastq data is fastqc. This tool will not 
perform any trimming or any other data processing. It will simply create a quality report in HTML format."""
software/fastqc sample1_r1.fastq sample1_r2.fastq sample1_out.R1.fq.gz sample1_out.R2.fq.gz

"""This is the most common way of visualising sequencing reads quality and most of you have been or will be
given a report like this."""
 
##################################         PART3      ###################################################       
###Alignment of the mitochndrial fastq sequence data to the mitochondrial reference sequence using BWA###
#########################################################################################################

##Alignment of the mitochndrial fastq sequence data to the mitochondrial reference sequence using BWA

"""BWA is a software package for mapping DNA sequences against a large reference genome, such as the 
human genome."""

"""Illumina/454/IonTorrent paired-end reads longer than ~70bp we use BWA MEM algorithm to align the fastq 
data to the reference genome"""

"""We will align our cleaned fastq files to the reference mitochondrial genome using BWA with the following
command;"""  

software/bwa mem reference/chrM.fa sample1_out.R1.fq.gz sample1_out.R2.fq.gz -R '@RG\tID:sample1\tSM:sample1\tLB:sample1\tPL:ILLUMINA' -o sample1.sam 

"""The software BWA mem has taken our cleaned fastq files, the indexed mitochondrial reference sequence
as arguments and generated the alignemnt file in sequence alignment format (.sam) as output"""

##converting the .sam output to the compressed common binary alignment format (.bam)

"""bam format has a lower data footprint than sam (smaller size on the disk), yet retains all of the same
information"""

software/samtools view -Sb sample1.sam -o sample1.bam

##sorting the .bam file

"""We sort the aligment file (.bam) by genomic coordinates so we will be able to then index it (next step).
Only coordinates-sorted alignment files can be indexed"""

software/samtools sort sample1.bam -o sample1_sorted.bam

##indexing the .bam file

"""Index a coordinate-sorted BAM or CRAM file for fast random access.""" 

software/samtools index sample1_sorted.bam

"""This command created a .bai file that is a companion file to the main .bam file. This file acts like an
external table of contents, and allows programs to jump directly to specific parts of the bam file without
reading through all of the sequences. This step will make the donwstream steps significantly faster."""

##deduplicate with picard tools

"""The deduplication step locates and tags duplicate reads in a BAM or SAM file, where duplicate reads are
defined as originating from a single fragment of DNA. Duplicates can arise during sample preparation e.g.
library construction using PCR."""

java -jar software/gatk-package-4.0.4.0-local.jar MarkDuplicates -I sample1_sorted.bam -O sample1_sorted_unique.bam \
	-M sample1_picard_metrics.txt
 
"""Duplicated reads will NOT be removed, they will just be tagged and ignored by downstream
steps."""

#base quality score recalibrations with GATK

"""We perform base quality score recalibration (BQSR) to recalculate the base calling scores
(it is a per sequenced base score) that the sequencer has determined. It's an important step as the short
variant calling algorithms used downstream rely heavily on the quality score assigned to the individual base
calls in each sequence read. BQSR comprises of two steps:"""

"""First step: BaseRecalibrator. This tool generates a recalibration table based on various covariates such as 
read group, reported quality score, machine cycle, and nucleotide context. The file (table) created on this step
will be used by the next step. Therefore, no changes will be made to our bam file on this step."""

java -jar software/gatk-package-4.0.4.0-local.jar BaseRecalibrator -I sample1_sorted_unique.bam -R reference/chrM.fa \
	--known-sites software/common_all_chrM.vcf.gz -O recal_data_table.txt
	
	
"""Second step: Apply Base Recalibration. This is the main and final step of BQSR.  The tool recalibrates the 
base qualities of the input reads based on the recalibration table produced on the previous step and outputs a 
new recalibrated BAM file."""	

java -jar software/gatk-package-4.0.4.0-local.jar ApplyBQSR -I sample1_sorted_unique.bam -R reference/chrM.fa \
	--bqsr-recal-file recal_data_table.txt -O sample1_sorted_unique_recalibrated.bam


##################################         PART4      ###################################################       
######## Identifying single nuclotide variants and small indels in our aligned mitochndrial data ########
#########################################################################################################

"""This is the part of the pipeline that we will use the alignment data to call variants, i.e., differences 
between the aligned reads and the used genome reference"""

##Variant Calling with GATK

"""GATK pipeline calls variants (SNPs and small INDELS simultaneously) using a tool called HaplotypeCaller.
This tool look through the alignments for regions with signs of variation (active regions). When it finds
an active region region like this, it discards the existing alignment and completely reassembles the reads
in that region making the calling more accurate."""

java -jar software/gatk-package-4.0.4.0-local.jar HaplotypeCaller -I sample1_sorted_unique_recalibrated.bam \
	-R reference/chrM.fa -G StandardAnnotation \
	-bamout sample1_sorted_unique_recalibrated.bamout.bam \
	-O sample1.vcf

"""Filter called variants. We will use GATK's VariantFiltration tool to filter the variants called previously.
Here we hard-filter variants based on certain criteria. In this particular example we will use:
-Genotype Quality less than 30.0: GQ < 30 to tag variants with low genotype quality
-Quality of Depth less than 1.5: QD < 1.5 to tag variants with low quality of depth
-Approximate read depth less than 6: DP < 6 to tag variants with low read coverage
-Strand Odds Ratio more than 10: SOR > 10 to tag variants with strand bias

This tool will TAG and NOT REMOVE the variants that false the Quality Control (fullfill at least one of
the above criteria)."""

java -jar software/gatk-package-4.0.4.0-local.jar VariantFiltration -V sample1.vcf -R reference/chrM.fa \
	--genotype-filter-expression "GQ < 30.0" --genotype-filter-name "LowGQ" \
	--filter-expression "QD < 1.5" --filter-name "LowQD" \
	--filter-expression "DP < 6" --filter-name "LowCoverage" \
	--filter-expression "SOR > 10.0" --filter-name "StrandBias" \
	-O sample1.filtered.vcf
	
"""The filtered/tagged variants will have filter-specific information on the FILTER column of the 
output vcf file. The variants that passed the filters will have a PASS label on their FILTER column.
Similarly to the deduplication step the concept here is to tag the problematic entries rather than 
to delete them."""
 
##################################         PART5      ###################################################       
######################### Functional annotation of the called genetic variants  #########################
#########################################################################################################

"""This is the part of the pipeline that we will functionally annotate the called genetic variants, an 
essential step to try predict the effect or function of an individual variant.

Given a list of variants and their genomic coordinates we will perform:

Gene-based annotations: Identify which genes and proteins will be affected by the reported variant 
	Example databases: RefSeq genes, UCSC genes, ENSEMBL genes, GENCODE genes, AceView genes, ... .

Region-based annotations: Identify specific genomic regions possibly affected by the reported variant
such as promoters, transcription f\actor binding sites, DNAse I hypersensitivity sites and others. 

Filter-based annotations: Cross listed variants with specific databases containing information about
them. For example, check which variants are reported in gnomAD populations and what's the allele  
frequency for these variants in this database or calculate a CADD score for each variant.
	Other example databases: 1000 Genome Project, Exome Aggregation Consortium (ExAC), 
	Genome Aggregation Database (gnomAD), CADD score, PolyPhen score. 

The goal is to create a single file containing several columns with annotation information for each 
variant.
"""

##################################         PART6      ###################################################       
############################## Aligment & Variants visualisation with IGV ###############################
#########################################################################################################

#add the code!!!







