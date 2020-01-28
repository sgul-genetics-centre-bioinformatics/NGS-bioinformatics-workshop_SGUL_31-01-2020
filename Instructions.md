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
