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
**Download the human mitochondrial refernce sequence:**
```bash
mkdir reference
wget --timestamping 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes/chrM.fa.gz' -O reference/chrM.fa.gz
```
**Unzip the compressed chrM fasta txt file**
```bash
gunzip reference/chrM.fa.gz
```
**View the top 10 lines of the file with the linux 'head' commands**
```bash
head reference/chrM.fa
```
**Prepare the FASTA file for use as a referenence for BWA**  
For the aligment algorithms of BWA, we first need to construct the FM-index for the reference genome 
(the index command):
```bash
software/bwa index reference/chrM.fa
```
**Prepare the FASTA file for use as a referenence for the [genome analysis toolkit (GATK)](https://gatk.broadinstitute.org/hc/en-us/articles/360036194592-Getting-started-with-GATK4)**  
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
