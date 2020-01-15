

##################################         PART1      ###################################################       
###Preparing our mitochondrial reference sequence fasta###
#########################################################################################################


##Download the human mitochondrial refernce sequence:
wget --timestamping 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg38/chromosomes/chrM.fa.gz' -O chrM.fa.gz

## unzip the compressed chrM fasta txt file
gunzip chrM.fa.gz

##view the top 10 lines of the file with the linux 'head' commands
head chrM.fa




###Prepare th FASTA file for use as a referenence for BWA###







###Prepare th FASTA file for use as a referenence for GATK###

"""Creat two needed files to access and safety check access to the reference files: a .dict dictionary of the contig names and sizes and a .fai fasta index file to allow efficient random access to the reference bases. 
You have to generate these files in order to be able to use a Fasta file as reference."""

##Creating the fasta index file

"""We use the faidx command in samtools to prepare the fasta index file. 
This file describes byte offsets in the fasta file for each contig, allowing us to compute exactly where a particular reference base at contig:pos is in the fasta file."""

samtools faidx chrM.fa

"""This produces a text file with one record per line for each of the fasta contigs. Each record is of the: contig, size, location, basesPerLine, bytesPerLine."""

##Creating the sequence dictionary file

"""We use CreateSequenceDictionary.jar from Picard to create a .dict file from a fasta file."""

java -jar CreateSequenceDictionary.jar R=chrM.fa O=chrM.dict

"""This produces a SAM-style header file describing the contents of our fasta file."""



##################################         PART2      ###################################################       
###Alignment of the mitochndrial fastq sequence data to the mitochondrial reference sequence using BWA###
#########################################################################################################

