# SGUL Workshop: Next Generation Sequencing data analysis
Friday, January 31, 2020, 14:00-17:00  
[St George's, University of London](https://www.sgul.ac.uk/)  
Room H5.2

## Organisers: 
- [Dr Alan Pittman](https://github.com/alanmichaelpittman100), Lecturer in Bioinformatics, SGUL
- [Dionysios Grigoriadis](https://github.com/digrigor), Bioinformatician, SGUL
- [SGUL Bioinformatics Unit](http://bioinformatics.sgul.ac.uk/)
- [SGUL Genetics Centre Bioinformatics](https://github.com/sgul-genetics-centre-bioinformatics)

## Next Generation Sequencing data analysis workshop
This hands-on beginners workshop, led by Dr Alan Pittman and Dionysios Grigoriadis, 
will cover the fundamental steps of analysing next-generation sequencing data; 
from processing, quality control and aligning raw sequence data to calling SNVs 
(short germline variants (Single Nucleotide Polymorphisms & short Indels) to obtain 
a reliable set of called variants for genetic analysis.  

Specifically:
-	Analysis of the human mitochondrial DNA (mtDNA) as an example.
-	Genome reference download and index.
-	Raw reads pre-processing and quality control.
-	Reads alignment to reference genome.
-	Alignment quality control and refinement.
-	Variant Calling and Filtration.
-	Annotation of the called variants.
-	Visualisation of alignments and called variants.

## Hands-on Tutorial
### Part 1 - Preparing our mitochondrial reference sequence fasta files
#### Download the human mitochondrial refernce sequence:
```bash
mkdir reference
wget --timestamping 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes/chrM.fa.gz' -O reference/chrM.fa.gz
```
---
#### Unzip the compressed chrM fasta txt file
```bash
gunzip reference/chrM.fa.gz
```
---
#### View the top 10 lines of the file with the linux 'head' commands
```bash
head reference/chrM.fa
```
---
#### Prepare the FASTA file for use as a referenence for BWA  
For the aligment algorithms of BWA, we first need to construct the FM-index for the reference genome 
(the index command):
```bash
software/bwa index reference/chrM.fa
```
---
#### Prepare the FASTA file for use as a referenence for the [genome analysis toolkit (GATK)](https://gatk.broadinstitute.org/hc/en-us/articles/360036194592-Getting-started-with-GATK4)  
Create two needed files to access and safety check access to the reference files:
- a .DICT dictionary of the contig names and sizes 
- a .FAI fasta index file to allow efficient random access to the reference bases.
You have to generate these files in order to be able to use a Fasta file as reference:

Creating the fasta index file. We use the faidx command in samtools to prepare the fasta index file. 
This file describes byte offsets in the fasta file for each contig, allowing us to compute exactly where 
a particular reference base at contig:pos is in the fasta file:

```bash
software/samtools faidx reference/chrM.fa
```
This produces a text file with one record per line for each of the fasta contigs. Each record is of the: 
contig, size, location, basesPerLine, bytesPerLine.

Creating the sequence dictionary file 
We use CreateSequenceDictionary.jar from Picard tools to create a .dict file from a fasta file:
```bash
java -jar software/gatk-package-4.0.4.0-local.jar CreateSequenceDictionary -R reference/chrM.fa -O reference/chrM.dict
```
This produces a SAM-style header file; simply describing the contents of our fasta file.

### Part 2 - Quality control of our mitochondrial NGS sequence data prior to alignments
#### View and inspect a FASTQ file
The pile of reads coming out of the sequencer is stored in FASTQ files. As a typical paired-end
sequencing experiment, each sequenced sample will produce two FASTQ files; One for the first pair of reads
and one for the second pair of reads.  

FASTQ files are usually quite big files containing information for all the reads produced. To get an insight
on how a FASTQ file look like we can simply open its first 12 lines:
```bash
head -n 12 sample1_r1.fastq
```
You can also try this command to scroll and navigate through the fastq file
```bash
more sample1_r1.fastq
```
---
#### Trimming by quality and adapter removal with fastp
This step can achieve 2 main aims:
With the program [fastp](https://github.com/OpenGene/fastp) we can: 
1. Assess the quality of our sequence fastq data, which we can visualise in 
the output HTML report file (which you can open with an internet browser).
2. Trim off low quality trailing bases from our sequence reads and detect 
and remove and remaining unwanted adapter sequencenes

To do that:
```bash
software/fastp -i sample1_r1.fastq -I sample1_r2.fastq -o sample1_out.R1.fq.gz -O sample1_out.R2.fq.gz --html \
 sample1_results.html --json sample1_results.json --report_title sample1_results
```

#### Alternative way of visualising the quality control of the reads with [FASTQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
An alternative way of visualising the quality of our sequence fastq data is fastqc. This tool will not 
perform any trimming or any other data processing. It will simply create a quality report in HTML format for the raw and the filtered reads.
```bash
software/fastqc sample1_r1.fastq sample1_r2.fastq sample1_out.R1.fq.gz sample1_out.R2.fq.gz
```
This is the most common way of visualising sequencing reads quality and most of you have been or will be
given a report like this. To learn more click [here](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

### Part 3 - Alignment of the mitochndrial fastq sequence data to the mitochondrial reference sequence using [BWA](http://bio-bwa.sourceforge.net/)
#### Alignment of the mitochndrial fastq sequence data to the mitochondrial reference sequence using BWA
[BWA](http://bio-bwa.sourceforge.net/) is a software package for mapping DNA sequences against a large reference genome, such as the 
human genome.  

Illumina/454/IonTorrent [paired-end reads](https://emea.illumina.com/science/technology/next-generation-sequencing/plan-experiments/paired-end-vs-single-read.html) longer than ~70bp we use BWA MEM algorithm to align the fastq 
data to the reference genome

We will align our cleaned FASTQ files to the reference mitochondrial genome using BWA with the following
command:
```bash
software/bwa mem reference/chrM.fa sample1_out.R1.fq.gz sample1_out.R2.fq.gz -R '@RG\tID:sample1\tSM:sample1\tLB:sample1\tPL:ILLUMINA' -o sample1.sam 
```
The software BWA mem has taken our cleaned fastq files, the indexed mitochondrial reference sequence
as arguments and generated the alignemnt file in sequence alignment format (.SAM) as output

---
#### Converting the .SAM output to the compressed common binary alignment format (.BAM)
BAM format has a lower data footprint than sam (smaller size on the disk), yet retains all of the same
information
```bash
software/samtools view -Sb sample1.sam -o sample1.bam
```
#### Sorting the .BAM file
We sort the aligment file (.BAM) by genomic coordinates so we will be able to then index it (next step).
Only coordinates-sorted alignment files can be indexed
```bash
software/samtools sort sample1.bam -o sample1_sorted.bam
```

#### Indexing the .BAM file
Index a coordinate-sorted BAM or CRAM file for fast random access.
```bash
software/samtools index sample1_sorted.bam
```
This command created a .bai file that is a companion file to the main .bam file. This file acts like an
external table of contents, and allows programs to jump directly to specific parts of the bam file without
reading through all of the sequences. This step will make the donwstream steps significantly faster.

---
#### Deduplication the aligned reads with picard tools
The deduplication step locates and tags duplicate reads in a BAM or SAM file, where duplicate reads are
defined as originating from a single fragment of DNA. Duplicates can arise during sample preparation e.g.
library construction using PCR. Read more [here](https://gatk.broadinstitute.org/hc/en-us/articles/360036350292-MarkDuplicates-Picard-)
```bash
java -jar software/gatk-package-4.0.4.0-local.jar MarkDuplicates -I sample1_sorted.bam -O sample1_sorted_unique.bam \
	-M sample1_picard_metrics.txt
```

---
#### Base quality score recalibrations with GATK
We perform base quality score recalibration (BQSR) to recalculate the base calling scores
(it is a per sequenced base score) that the sequencer has determined. It's an important step as the short
variant calling algorithms used downstream rely heavily on the quality score assigned to the individual base
calls in each sequence read. BQSR comprises of two steps:
1. BaseRecalibrator. This tool generates a recalibration table based on various factors of the aligned reads such as 
read group, reported quality score, machine cycle, and nucleotide context. The file (table) created on this step
will be used by the next step. Therefore, no changes will be made to our bam file on this step.
```bash
java -jar software/gatk-package-4.0.4.0-local.jar BaseRecalibrator -I sample1_sorted_unique.bam -R reference/chrM.fa \
	--known-sites software/common_all_chrM.vcf.gz -O recal_data_table.txt
```
1. Apply Base-Recalibration. This is the main and final step of BQSR.  The tool recalibrates the 
base qualities of the input reads based on the recalibration table produced on the previous step and outputs a 
new recalibrated BAM file:
```bash
java -jar software/gatk-package-4.0.4.0-local.jar ApplyBQSR -I sample1_sorted_unique.bam -R reference/chrM.fa \
	--bqsr-recal-file recal_data_table.txt -O sample1_sorted_unique_recalibrated.bam
```


