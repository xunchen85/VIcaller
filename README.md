# VIcaller
A software to detect virome-wide viral integrations

## 1 Introduction
Viral Integration caller (VIcaller) is a bioinformatics tool designed for identifying viral integration events using high-throughput sequencing (HTS) data. VIcaller is developed under Linux platform. It uses both FASTQ files or aligned BAM files as input. It also supports both single-end and paired-end reads. VIcaller contains one main Perl script, VIcaller.pl, that include three main functions: 1) detect, which will detect virome-wide candidate viruses and integration events; 2) validate, which will perform the in silico validation on those candidate viral integrations; 3) calculate, which will calculate the integration allele fraction. We also generated a comprehensive viral reference genome library with 411,195 unique whole and partial genomes, covering all six virus taxonomic classes. The viral reference genome library  also comes with a taxonomy database in a defined format that give the virus name, etc. 

## 2 Availability
VIcaller is an open-source software. VIcaller.v1.1 source code is available at www.uvm.edu/genomics/software/VIcaller, virome-wide library is also available here, and the human reference genome (GRCh38) can be obtained from UCSC website.

## 3 VIcaller installation
### 3.1 Unzip the VIcaller installer
#### Unzip the installer and change the directory
```
$ tar vxzf VIcaller.tar.gz
$ cd VIcaller/
$ mkdir Tools
```
#### 3.2 Install the dependent Perl libraries and tools
a) Currently VIcaller relies on the following dependencies to be compiled (contact Dr. Xun Chen if you need help get those tools or Perl libraries installed).  
b) Obtain the installed file from the following links.  
c) Follow the instruction to successfully install each tool (contact server manager if there is any compile issues).  
d) Check or install the listed Perl libraries using cpan, cpanm or other methods.  

#### Install each of the listed tools
•	BWA (default version: v0.7.10): https://github.com/lh3/bwa/tree/master/bwakit  
•	Bowtie2 (default version: v2.2.7): https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.7/  
•	TopHat2 (v2.1.1): http://ccb.jhu.edu/software/tophat/index.shtml  
•	BLAT (default version: v.35): http://genomic-identity.wikidot.com/install-blat  
•	BLAST+ (default version: v2.2.30): http://mirrors.vbi.vt.edu/mirrors/ftp.ncbi.nih.gov/blast/executables/blast%2B/2.2.30/  
•	SAMtools (default version: v1.6): https://sourceforge.net/projects/samtools/  
•	HYDRA (default version: 0.5.3): https://code.google.com/archive/p/hydra-sv/downloads  
•	NGS QC Toolkit (default version: v2.3.3): http://genomic-identity.wikidot.com/install-blat  
a) Copy the script “TrimmingReads_sanger.pl” under the VIcaller/Scripts/ folder to the installed NGSQCToolkit_v2.3.3/Trimming/ folder  
•	FastUniq (Default version: v1.1): https://sourceforge.net/projects/fastuniq/  
•	SE-MEI (modified): https://github.com/dpryan79/SE-MEI (original version), the modified version can be found under the VIcaller/Scripts/ folder  
a) Copy the modified SE-MEI installer (SE-MEI-master.tar.gz) under the VIcaller/Scripts/ folder to the VIcaller/Tools/ folder  
b) Install the modified SE-MEI tool follow the README file  
•	RepeatMasker (default version: v4.0.5):  
a) Install RepeatMasker: http://www.repeatmasker.org/  
b) Install RMBlast aligner: http://www.repeatmasker.org/RMBlast.html  
c) Compile the Repbase database: https://www.girinst.org/repbase/   
•	MEME (default version: v4.11.1): http://web.mit.edu/meme_v4.11.4/share/doc/download.html  
•	TRF (default version: v4.07b): https://tandem.bu.edu/trf/trf.html  

#### Install Perl libraries
```
$ cpan String::Approx
$ cpan Time::HiRes
$ cpan Test::Most
$ cpan Bio::Seq
$ cpan Bio::SeqIO
$ cpan Bio::DB::GenBank
$ cpan IO::Zlib
```

###	3.3 Prepare databases
#### Index human reference genome using BWA, Bowtie2, and BLAST+ separately:
```
$ bwa index -a bwtsw hg38.fa
$ bowtie2-build hg38.fa hg38.fa
$ makeblastdb -in hg38.fa -dbtype nucl
```

#### Index viral database using BWA, Bowtie2, and BLAST+ separately:
```
$ bwa index -a bwtsw hg38.fa
$ bowtie2-build hg38.fa hg38.fa
$ makeblastdb -in hg38.fa -dbtype nucl
```

### 3.4 Prepare the VIcaller config file
#### Example of VIcaller.config
```
export PERL5LIB=/users/xchen/.cpan/build/
export PATH=$PATH:/users/xchen/VIcaller/Tools/bowtie2-2.2.7/
# human_genome = /users/xchen/VIcaller/Database/Human/human_genome_hg38/ hg38.fa
# human_genome_tophat = /users/xchen/VIcaller/Database/Human/human_genome_hg38/ hg38.fa
# virus_genome = /users/xchen/VIcaller/Database/Virus/virus_db_090217/virus_db_090217.fa          
# virus_taxonomy = /users/xchen/VIcaller/Database/Virus/virus_db_090217/virus_db_090217.taxonomy
# virus_list = /users/xchen/VIcaller/Database/Virus/virus_db_090217/virus_db_090217.virus_list
# vector_db = /gpfs2/dli5lab/CAVirus/Database/Vector/Vector.fa
# cell_line = /users/xchen/VIcaller/Database/cell_line.list
# bowtie_d = /users/xchen/VIcaller/Tools/bowtie2-2.2.7/
# tophat_d = /users/xchen/VIcaller/Tools/tophat-2.1.1.Linux_x86_64/
# bwa_d = /users/xchen/VIcaller/Tools/bwa-master/
# samtools_d = /users/xchen/VIcaller/Tools/samtools-1.6/
# repeatmasker_d = /users/xchen/VIcaller/Tools/RepeatMasker/
# meme_d = /users/xchen/VIcaller/Tools/meme_4.11.1/
# NGSQCToolkit_d = /users/xchen/VIcaller/Tools/NGSQCToolkit_v2.3.3/
# fastuniq_d = /users/xchen/VIcaller/Tools/FastUniq/
# SE_MEI_d = /users/xchen/VIcaller/Tools/SE-MEI/
# hydra_d = /users/xchen/VIcaller/Tools/Hydra-Version-0.5.3/
# blat_d = /users/xchen/bin/x86_64/
# blastn_d = /users/xchen/VIcaller/Tools/ncbi-blast-2.2.30+-src/
```

#### Check the generated VIcaller.config file
#. Make sure the space between “#” and parameters.  
#. Make sure the directory for the Perl library is correct or the libraries are available in the path if you install them locally.  
#. Make sure the Bowtie2 directory is correct or it is available in the path (recommended) if you are going to analyze RNA-seq data.  
#. Make sure the human and virus databases existed and correctly indexed.  

## 4 VIcaller command line
```
$ perl VIcaller.pl <functions> [arguments]
```

### 4.1 Detect candidate viral integrations
#### Command line
```
$ perl VIcaller.pl detect [arguments]
```
#### Examples
a) WGS data in single-end fastq format:  
```
$ perl VIcaller.pl detect -d WGS -i seq -f .fastq.gz -s single-end -t 12
```
b) RNA data in paired-end fastq format (set bowtie2 path before run the following command):
```
$ perl VIcaller.pl detect -d RNA-seq -i seq -f .fastq.gz -s paired-end -t 12
```
c) RNA alignment data in bam format  (Note: Human reference genome should be the same as the bam file)
```
$ perl VIcaller.pl detect -d RNA-seq -i seq -f .bam -s paired-end -t 12
```
### 4.2 Validate candidate viral integrations
#### Command line
```
$ perl VIcaller.pl validate [arguments]
```
#### Example
```
$ perl VIcaller.pl validate -i seq -S seq_1_24020575_24020787_HPV16_218931404 -G 218931404 -V HPV16
```

### 4.3 Calculate allele fraction
#### Command line
```
$ perl VIcaller.pl calculate [arguments]
```

#### Example
```
$ perl VIcaller.pl calculate -i seq -f .fastq.gz -S -C 1 -P 24020575 -B 2 -N 20
```

## 5 Output
#### Output and file list
The candidate viral integrations detected by VIcaller are kept in the file with suffix of “.output” in Viral integration Format (VIF), with the visualization of the aligned read sequences in the file with suffix of “.visualization”. After in silico validation and allele fraction calculation, the results are also kept in the output file. “seq” is an example sample ID.

#### Header of the output file
Details are included in the Manual PDF file

## 6 FAQ
### 6.1 How to annotate the detected viral integrations?
        The following Linux command can be used to extract the information required to run human genome functional annotation tools. The VIcaller output file is “seq.output”, and for example, if the functional annotation software is SnpEff, the following command line will extract the information required to run SnpEff. The output from using this command will be the input file for SnpEff. 
```
$ awk '{if ($7!="Chr.")print$7"\t"$17"\t.\tA\tT\t."}' seq.output
```
### 6.2 What is the difference between “Fast” mode and “Standard” mode?
       “Fast” mode is significantly faster than “Standard” mode. However, the “Fast” mode does not analyze viral reads, which are supporting evidence for distinguishing between viral integrations and viral infections.  

### 6.3 How to use the viral integration data from VIcaller for integration enrichment analysis?
        VIcaller analyzes individual samples and then generates a list of viral integrations for each sample. Viral integration enrichment (bias) analysis, which is a statistical analysis, requires inclusion of a group of samples. The enrichment analysis has to be performed solely. There are multiple statistical models to calculate/determine enrichment hotspots (such as simulation-based Z score test). There are many available tools and R packages that can be selected for enrichment analysis. Users may have different preference on statistical models to fit in their actual samples/data. 

### 6.4 Can I use the published tools that were designed for detecting transposable element insertions to identify virome-wide integrations?
        VIcaller uses the reads that are commonly used in transposable element insertion and other structural variation detection tools. However, because VIcaller is specifically designed to identify virome-wide integrations, it has significant advantages over alignment-based transposable element insertion detection tools for viral integration analysis, which are designed to extract and mainly use (human’s) anomalous reads specifically. For example, 1) VIcaller supports the use of virome-wide library as the reference to detect any characterized viruses, while most transposable element detection tools use transposable element sequences as the reference; and 2) VIcaller implements viral integration-specific quality control procedures and implements additional steps to in silico verify detected viral integrations. We have tried to compare VIcaller with other transposable element insertion detection software, e.g., MELT. MELT failed to run in a virome-wide fashion after we replaced MELT’s default consensus transposable element reference sequences with our virome-wide database. We further tested whether MELT was able to detect simulated candidate viral integrations, and we found that although MELT did run, it was not able to detect any of these integrations.

### 6.5 Where to get the human reference genome? 
        You can download hg38 here: http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/. It is recommended to use the latest hg38.fa.gz file for indexing.

### 6.6 Can I use other virome-wide libraries?
        You can use other viral databases as the reference. However, the final output may not include the viral names or other taxonomy information. The reads that multiple-mapped different viral sequences from the same virus may not be efficiently recovered for the detection of viral integrations.

## Citation
Xun Chen, Jason Kost, Arvis Sulovari, Nathalie Wong, Winnie S. Liang, Jian Cao, and Dawei Li. A virome-wide clonal integration analysis platform for discovering cancer viral etiologies. Under review.

## Download
www.uvm.edu/genomics/software/VIcaller

## Copyright
VIcaller is licensed under the Creative Commons Attribution-NonCommercial 4.0 International license. It may be used for non-commercial use only. For inquiries about a commercial license, please contact the corresponding author at dawei.li@uvm.edu or The University of Vermont Innovations at innovate@uvm.edu.


