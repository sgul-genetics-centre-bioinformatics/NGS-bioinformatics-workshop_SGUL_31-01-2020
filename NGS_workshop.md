![alt text](https://ukeducationguide.com/wp-content/uploads/2014/10/stgeorgeslondon.jpg "St George's, University of London") 
# SGUL Workshop: Next Generation Sequencing data analysis
Friday, January 31, 2020, 14:00-17:00  
[St George's, University of London](https://www.sgul.ac.uk/)  
Room H5.2

## Organisers: 
- [Dr Alan Pittman](https://github.com/alanmichaelpittman100), Lecturer in Bioinformatics, SGUL  
	[✉ apittman@sgul.ac.uk](mailto:apittman@sgul.ac.uk?subject=SGUL%2Workshop)
- [Dionysios Grigoriadis](https://github.com/digrigor), Bioinformatician, SGUL  
	[✉ dgrigori@sgul.ac.uk](mailto:dgrigori@sgul.ac.uk?subject=SGUL%2Workshop)
- [SGUL Bioinformatics Unit](http://bioinformatics.sgul.ac.uk/)
- [SGUL Genetics Centre Bioinformatics](https://github.com/sgul-genetics-centre-bioinformatics)

## Workshop set-up
For this workshop we are going to use the command line enviroment of the University's STATS3 server, which is accessible to all staff/students with valid SGUL credentials.

Instructions:
1. Open File Explorer, locate the N: drive, find and run MobaXterm_Personal_10.4.exe
2. Click "Start Local Terminal"
3. Login by typing:
	```bash
	ssh yourusername@stats3.sgul.ac.uk
	```
	Type your password when prompted and press Enter (Note that your password will not be printied to the screen as you type).

4.    
	**IF YOU ARE PHD STUDENT OR STAFF**:    
	This login defaults to your home (H:) drive where we will be working.  
	Move on to next step.

	**IF YOU ARE UNDERGRADUATE OR POSTGRADUATE STUDENT**:  
	This login defaults to your home (H:) drive where we will be working.  
	Move inside your working directory by typing:
	```bash
	/homedirs8/workshop/<username>
	```

5. Type:
	```bash
	git clone https://github.com/sgul-genetics-centre-bioinformatics/NGS-bioinformatics-workshop_SGUL_31-01-2020.git
	```

6. Wait for the directory to download and then move inside the directory:
	```bash
	cd NGS-bioinformatics-workshop_SGUL_31-01-2020
	```

7. You can view a list of all the contents of the directory you are in, by typing:
	```bash
	ls
	```

You are now ready to start this workshop! Relax and enjoy :)

## Learning Objectives
- Quick overview of next generation sequencing 
- Linux command line
- Next-generation sequencing software/tools 
	- Genome reference download and index
	- Raw sequencing reads pre-processing and quality control
	- Reads alignment to the reference genome
	- Alignment quality control and refinement
	- Variant calling and filtering
	- Annotation of the called variants
	- Visualisation of alignments and called variants

## Hands-on Tutorial

---
### Part 1 - Preparing our mitochondrial reference sequence fasta files
#### Download the human mitochondrial reference sequence:
```bash
mkdir reference
wget --timestamping 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes/chrM.fa.gz' -O reference/chrM.fa.gz
```
---
#### Unzip the compressed chrM fasta text file
```bash
gunzip reference/chrM.fa.gz
```
---
#### View the top 10 lines of the file with the linux 'head' command
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
#### Prepare the FASTA file for use as a reference for the [genome analysis toolkit (GATK)](https://gatk.broadinstitute.org/hc/en-us/articles/360036194592-Getting-started-with-GATK4)  
Create two files needed to access and safety check access to the reference files:
- a .DICT dictionary of the contig names and sizes 
- a .FAI fasta index file to allow efficient random access to the reference bases.
You have to generate these files in order to be able to use a Fasta file as reference:

Creating the fasta index file. We use the faidx command in samtools to prepare the fasta index file. 
This file describes byte offsets in the fasta file for each contig, allowing us to compute exactly where 
a particular reference base at contig:pos is in the fasta file:

```bash
software/samtools faidx reference/chrM.fa
```
This produces a text file with one record per line for each of the fasta contigs. Each record is of the format: 
contig, size, location, basesPerLine, bytesPerLine.

Creating the sequence dictionary file. We use CreateSequenceDictionary.jar from Picard tools to create a .dict file from a fasta file:
```bash
java -jar software/gatk-package-4.0.4.0-local.jar CreateSequenceDictionary -R reference/chrM.fa -O reference/chrM.dict
```
This produces a SAM-style header file; simply describing the contents of our fasta file.

---
### Part 2 - Quality control of our mitochondrial NGS sequence data prior to alignments
#### View and inspect a FASTQ file
The pile of reads coming out of the sequencer is stored in FASTQ files. As a typical paired-end
sequencing experiment, each sequenced sample will produce two FASTQ files; One for the first in the pair 
and one for the second in each pair of reads.  

FASTQ files are usually quite big files containing information for all the reads produced. To get an insight
into how a FASTQ file looks, we can simply open the first 12 lines:
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

To view the generated .HTML report with firefox type the command below and wait for a few seconds.
```bash
firefox sample1_results.html
```

If this command doesn't work, you can navigate to the sample1_results.html file through your Windows environment and open the file with a browser (Google Chrome, Mozilla Firefox, etc...)
#### Alternative way of visualising the quality control of the reads with [FASTQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
An alternative way of visualising the quality of our sequence fastq data is with fastqc. This tool will not 
perform any trimming or any other data processing. It will simply create a quality report in HTML format for the raw and the filtered reads.
```bash
software/fastqc sample1_r1.fastq sample1_r2.fastq sample1_out.R1.fq.gz sample1_out.R2.fq.gz
```

To view the generated .HTML report with firefox type the command below and wait for a few seconds.
```bash
firefox sample1_out.R1_fastqc.html
```
If this command doesn't work, you can navigate to the sample1_out.R1_fastqc.html file through your Windows environment and open the file with a browser (Google Chrome, Mozilla Firefox, etc...)

This is the most common way of visualising sequencing reads quality and most of you have been or will be
given a report like this. To learn more click [here](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/).  

---
### Part 3 - Alignment of the mitochndrial fastq sequence data to the mitochondrial reference sequence using [BWA](http://bio-bwa.sourceforge.net/)
#### Alignment of the mitochndrial fastq sequence data to the mitochondrial reference sequence using BWA
[BWA](http://bio-bwa.sourceforge.net/) is a software package for mapping DNA sequences against a large reference genome, such as the 
human genome.  

For Illumina [paired-end reads](https://emea.illumina.com/science/technology/next-generation-sequencing/plan-experiments/paired-end-vs-single-read.html) longer than ~70bp we use the BWA MEM algorithm to align the fastq 
data to the reference genome

We will align our cleaned FASTQ files to the reference mitochondrial genome using BWA with the following
command:
```bash
software/bwa mem reference/chrM.fa sample1_out.R1.fq.gz sample1_out.R2.fq.gz \
-R '@RG\tID:sample1\tSM:sample1\tLB:sample1\tPL:ILLUMINA' -o sample1.sam 
```
The software BWA mem has taken our cleaned fastq files, the indexed mitochondrial reference sequence
as arguments and generated the alignment file in sequence alignment format (.SAM) as output

---
#### Converting the .SAM output to the compressed common binary alignment format (.BAM)
BAM format has a lower data footprint than SAM (smaller size on the disk), yet retains all of the same
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
2. Apply Base-Recalibration. This is the main and final step of BQSR.  The tool recalibrates the 
base qualities of the input reads based on the recalibration table produced on the previous step and outputs a 
new recalibrated BAM file:
```bash
java -jar software/gatk-package-4.0.4.0-local.jar ApplyBQSR -I sample1_sorted_unique.bam -R reference/chrM.fa \
--bqsr-recal-file recal_data_table.txt -O sample1_sorted_unique_recalibrated.bam
```
---
### Part 4: Variant Calling: Identifying single nucleotide variants and small indels in our aligned mitochndrial data
This is the part of the pipeline that we will use the alignment data to call variants, i.e., differences 
between the aligned reads and the reference genome.

2 steps are essential for variant calling:
#### 1. HaplotypeCaller step
The GATK pipeline calls variants (SNPs and small INDELS simultaneously) using a tool called HaplotypeCaller.
This tool look through the alignments for regions with signs of variation (active regions). When it finds
an active region, it discards the existing alignment and completely reassembles the reads
in that region making the calling more accurate.
```bash
java -jar software/gatk-package-4.0.4.0-local.jar HaplotypeCaller -I sample1_sorted_unique_recalibrated.bam \
-R reference/chrM.fa -G StandardAnnotation \
-bamout sample1_sorted_unique_recalibrated.bamout.bam \
-O sample1.vcf
```

#### 2. Filtering of the called variants
Filter called variants. We will use GATK's VariantFiltration tool to filter the variants called previously.
Here we hard-filter variants based on certain criteria. In this particular example we will use:
-Genotype Quality less than 30.0: GQ < 30 to tag variants with low genotype quality
-Quality of Depth less than 1.5: QD < 1.5 to tag variants with low quality of depth
-Approximate read depth less than 6: DP < 6 to tag variants with low read coverage
-Strand Odds Ratio more than 10: SOR > 10 to tag variants with strand bias

This tool will TAG and NOT REMOVE the variants that fail the Quality Control (fullfill at least one of
the above criteria).
```bash
java -jar software/gatk-package-4.0.4.0-local.jar VariantFiltration -V sample1.vcf -R reference/chrM.fa \
--genotype-filter-expression "GQ < 30.0" --genotype-filter-name "LowGQ" \
--filter-expression "QD < 1.5" --filter-name "LowQD" \
--filter-expression "DP < 6" --filter-name "LowCoverage" \
--filter-expression "SOR > 10.0" --filter-name "StrandBias" \
-O sample1.filtered.vcf
```

The filtered/tagged variants will have filter-specific information in the FILTER column of the 
output vcf file. The variants that passed the filters will have a PASS label in their FILTER column.
Similarly to the deduplication step the concept here is to tag the problematic entries rather than 
to delete them.

---
#### Inspect the output vcf file sample1.filtered.vcf
To view the contents of the output filtered vcf file you can navigate through it.  
Apart from the head and more commands we saw earlier, you can also use the `less` command in
Linux to inspect the contents of a file:
```bash
less sample1.filtered.vcf
```
Or you can open it in the notepad editor from the Windows environment.

---
### Part 5 - Functional annotation of the called genetic variants
This is the part of the pipeline that we will functionally annotate the called genetic variants, an 
essential step to try to predict the effect or function of an individual variant.

Given a list of variants and their genomic coordinates (like the vcf file you just inspected) we will perform functional annotations 
for each one of them.

Variants are annotated with annotation databases. An annotation database can be one of the following types:
- Gene-based annotations: Identify which genes and proteins will be affected by the reported variant   
Example databases: RefSeq genes, UCSC genes, ENSEMBL genes, GENCODE genes, AceView genes, ... .

- Region-based annotations: Identify specific genomic regions possibly affected by the reported variant
such as promoters, transcription factor binding sites, DNAse I hypersensitivity sites and others. 

- Filter-based annotations: Cross listed variants with specific databases containing information about
them.  
For example, check which variants are reported in gnomAD populations and what's the allele  
frequency for these variants in this database or calculate a CADD score for each variant.
	Other example databases: 1000 Genome Project, Exome Aggregation Consortium (ExAC), 
	Genome Aggregation Database (gnomAD), CADD score, PolyPhen score.

We are going to use databases from all of the above three catergories to annotate our variants.

The goal is to create a single file containing several columns with annotation information for each 
variant.

#### Annotate online using the [Variant Effect Predictor (VEP)](https://www.ensembl.org/vep) by Ensembl
For this part we are not going to use a command line tool. We are going to use a web-based tool instead
which is called [Variant Effect Predictor (VEP)](https://www.ensembl.org/vep) by Ensembl

VEP also has a [command line stand alone version](http://www.ensembl.org/info/docs/tools/vep/script/vep_download.html) which can be freely downloaded and used. 

**Please follow the instructions below:**
1. Open an internet browser (Mozilla Firefox, Google Chrome).

2. Navigate to https://www.ensembl.org/Multi/Tools/VEP.

3. **! Double check the human genome build you are using ! (Here, we are using hg38/GRCh38)**.

4. Fill the job form using:
	- Name for this job: `sgul_workshop_` followed by your initials (e.g. `sgul_workshop_DG`).
	- Input data -> Or upload file: -> Choose file -> Navigate and select your `sample1.filtered.vcf` file.
	- Transcript database to use: Ensembl/GENCODE transcripts.
	- Make sure you will leave the Additional configurations field untouched using the pre-defined default fields.

4. Click RUN.

5. Wait for the job to finish.

6. When finished, click [View results].

---
### Part 6: Aligment & Variants visualisation with [Integrative Genomics Viewer (IGV)](https://igv.org/)

This is the part of the pipeline that we will visualise the alignments generated in **Part 3** and the filtered variant calls generated in **Part 4**.

For this part we are going to use the linux version of the [IGV application](https://igv.org/)

1. At the linux terminal type:
```bash
software/igv.sh
```
2. Open the IGV software (A window will normally open in your windows environment, please check the task bar).
3. Go to `Genomes` and then select `Load genomes from file...`
	* Locate the chrM.fa in your references directory and open.
4. Go to `File` and then select `Load from file...`
	* Locate the `sample1_sorted_unique_recalibrated.bam` and open.
  
You can also load a track showing the actual called variants:
- Go to `File` and then select `Load from file...`
- Locate the `sample1.filtered.vcf` and open.
  
  
  
There is also an [online version of IGV](https://igv.org/app/), which you can use without any installation.


# Quiz
If you didn't have enough already, time to spend some time by playing around with this Quiz!

Have fun and **Don't be shy**. It's normal to get stuck in many of these questions, so feel free to grab a demonstrator and seek advice! :)

---
1. In the VEP annotation results table, find the genomic coordinates of the missense variant with the highest Polyphen score.

---
2. Download the annotated vcf file using the VCF download button on the upper-right corner.

---
3. Get the genomic position of the variant you selected in question 1 and extract this variant’s information from the sample1.filtered.vcf file (using the [`grep` command](https://linuxtechlab.com/learning-grep-command-with-examples/)) and save it into a new file ([using the > symbol](https://askubuntu.com/questions/420981/how-do-i-save-terminal-output-to-a-file)).

---
4. By looking into the new file you created:
	- Find out if the sequenced sample is homozygous or heterozygous for this variant.
	- Report the Strand Odds Ratio for this variant.
	- Report the Genotype Quality and the Quality of Depth values for this variant.  
Hint! Check the .VCF specification [here](https://samtools.github.io/hts-specs/VCFv4.2.pdf) and Part4 of the workshop.

---
5. Load the IGV software with the alignment (.BAM) file and the variant calling (.VCF) file.
	- Point the IGV viewer to the previously selected variant.
	- Does it seem real to you? What do you think?  
Hint: Part 6 of the workshop.

---
6. Count the reads covering the position of the selected variants in the alignment file (.BAM) using the samtools tool.  
Hint! You can find samtools in software/samtools).  
Hint! Check samtools documentation here: http://www.htslib.org/doc/samtools.html.  
Hint! Use the samtools view function.  


