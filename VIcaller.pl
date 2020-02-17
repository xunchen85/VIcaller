#!usr/bin/env perl
#
# Author: 	Xun Chen
# Email: 	Xun.Chen@uvm.edu or xunchen85@gmail.com
# Date:		01/08/19
# Version:      v1.1
#
#
############# updates
#
# Version v1.1:
#	04/15/2019	Debug the validate function;
#	01/10/2019	Debug VIcaller using the RNA-seq FASTQ file as input;
# 	01/08/2019	Remove the function of installing depedent tools due to the links to get the installers expired;
#	12/02/2018	Debug VIcaller for the use of split reads;
# 	06/17/2018	Correct the input data format, now accept input files with any names;
#
# Version v1.0:
# 	07/22/2018	Adding the functions of annotating vector sequences;
# 	07/22/2018	Adding the functions of annotating cell line contaminations;
#
#
# Questions for helps:
#
# Contact Xun Chen Ph.D. for any questions

use Getopt::Long qw(:config no_ignore_case);
use strict;
use warnings;
use Cwd qw();

#### Parameter variables;
# main
my $function = $ARGV[0];
my $input_sampleID;
my $threads = 16;
my $help;
my $config;

# detect
my $data_type = "WGS";
my $mode = "standard";
my $file_suffix=".fq.gz";
my $sequencing_type="paired-end";
my $length_insertsize=1500;
my $directory=$0;
my $QS_cutoff;
my $repeat;
my $align_back_to_human;
my $build="hg38";

# validate
my $String;
my $GI;
my $Virus;

# calculate
my $File_suffix_bam;
my $Index_sort;
my $Chr;
my $Position;
my $Breakpoint;
my $Number_reads;


##### variables
GetOptions(
           "i|input_sampleID=s" => \$input_sampleID,
           "h|help" => \$help,
           "t|threads=i" => \$threads,
	   "c|config=s" => \$config,
           "f|file_suffix=s" => \$file_suffix,
           "d|data_type=s" => \$data_type, 
           "m|mode=s" => \$mode,
           "s|sequencing_type=s" => \$sequencing_type,
           "q|QS_cutoff=i" =>\$QS_cutoff,
	   "r|repeat" => \$repeat,
	   "a|align_back_to_human" => \$align_back_to_human,
           "l|length_insertsize=i" =>\$length_insertsize,
           "b|build=s" =>\$build,

           "S|String=s" => \$String,
           "G|GI=s" => \$GI,
           "V|Virus=s" => \$Virus,

           "F|File_suffix_bam=s" => \$File_suffix_bam,
	   "I|Index_sort" => \$Index_sort,
	   "C|Chr=s" => \$Chr,
           "P|Position=i" => \$Position,
           "B|Breakpoint=i" => \$Breakpoint,
           "N|Number_reads=i" => \$Number_reads,
           );

##### time
my @months = qw( Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec );
my @days = qw(Sun Mon Tue Wed Thu Fri Sat Sun);
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime();

print "\n\n\n#############################################################\n";
print "###############                              \n";
print "###############      $hour:$min:$sec, $months[$mon].$mday\n";
print "###############      Author: Xun Chen (Xun.Chen\@uvm.edu and xunchen85\@gmail.com)\n";
print "###############                              \n";
print "#############################################################\n\n\n\n";
print "~~~~~ Begin...\n";

my $input_sampleID2 =$String;
my $candidate_virus=$Virus;
my @groups;
my $ref;

##### default value
if (!defined ($QS_cutoff)) {
  $QS_cutoff = 20;
}


##### check if the function is defined
if(!($function)){
  prtErr("#Error: please defined the function\n");
  prtUsa();
  exit;
}
elsif($function ne "detect" && $function ne "calculate" && $function ne "config" && $function ne "validate" & !defined($help)) {
  prtErr("#Error: unknown function defined\n");
  prtUsa();
  exit;
} else {
  print "\n# Functions: VIcaller.pl $function..........\n\n";
}

##### check if sampleID provided
if(!defined($input_sampleID) & !defined($help) && $function eq "detect") {
  prtErr("#Error: no samples are provided");
  prtUsa();
  exit;
}

##### show help page
if(defined($help)){
  prtUsa();
  exit;
}

##### extract the path of the installed VIcaller software
my @directory=split /\//, $directory;
if($#directory>0){
  $directory=join("/",@directory[0..$#directory-1]);
  $directory=$directory."/";
} else {
  $directory="";
}
 ##### load_config();
 my $human_genome;
 my $human_genome_tophat;
 my $virus_genome;
 my $virus_taxonomy;
 my $virus_list;
 my $vector_db;
 my $cell_line;
 my $bowtie2_d="";
 my $tophat_d="";
 my $bwa_d="";
 my $samtools_d="";
 my $repeatmasker_d="";
 my $dust_d="";
 my $NGSQCToolkit_d="";
 my $fastuniq_d="";
 my $SE_MEI_d="";
 my $hydra_d="";
 my $blat_d="";
 my $blastn_d="";
 if (defined($config)){
   open CINFIG,"$config";
 } elsif (-e $directory."VIcaller.config") {
   open CINFIG,"${directory}VIcaller.config";
 } else {
   print "# Error: missing config file, you can either specify one, or create one using function \"config\"\n";
   exit;}
 while(<CINFIG>){
   my @line=split;
   if(!(@line)){next;}    
   if($line[0] eq "\#"){
     if($line[1] eq "human_genome"){$human_genome=$line[3];}
     elsif($line[1] eq "human_genome_tophat"){$human_genome_tophat=$line[3];}
     elsif($line[1] eq "virus_genome"){$virus_genome=$line[3];}
     elsif($line[1] eq "virus_taxonomy"){$virus_taxonomy=$line[3];}
     elsif($line[1] eq "virus_list"){$virus_list=$line[3];}
     elsif($line[1] eq "cell_line"){$cell_line=$line[3];}
     elsif($line[1] eq "vector_db"){$vector_db=$line[3];}
     elsif($line[1] eq "bowtie2_d"){$bowtie2_d=$line[3];}
     elsif($line[1] eq "tophat_d"){$tophat_d=$line[3];}
     elsif($line[1] eq "bwa_d"){$bwa_d=$line[3];}
     elsif($line[1] eq "samtools_d"){$samtools_d=$line[3];}
     elsif($line[1] eq "repeatmasker_d"){$repeatmasker_d=$line[3];}
     elsif($line[1] eq "dust_d"){$dust_d=$line[3];}
     elsif($line[1] eq "NGSQCToolkit_d"){$NGSQCToolkit_d=$line[3];}
     elsif($line[1] eq "fastuniq_d"){$fastuniq_d=$line[3];}
     elsif($line[1] eq "SE_MEI_d"){$SE_MEI_d=$line[3];}
     elsif($line[1] eq "hydra_d"){$hydra_d=$line[3];}
     elsif($line[1] eq "blat_d"){$blat_d=$line[3];}
     elsif($line[1] eq "blastn_d"){$blastn_d=$line[3];}
   }
 }
 system ("chmod a+x ${directory}VIcaller.config");
 system ("${directory}./VIcaller.config");

my $thread_1=$threads-1;
my $order=0;
#####  1. detect function
if ($function eq "detect"){
 check_input();
 print "~~~~~ step 1\n";
 ##### extract supporting reads
 if(${file_suffix} =~ "bam"){
   $order=1;
   convert_bamtofastq(${input_sampleID});
   align_to_hg(${input_sampleID}."_h1",".1fq");
   $order=2;
   convert_bamtofastq(${input_sampleID}."_h1");
   system("mv ${input_sampleID}_h1_sm.bam ${input_sampleID}_sm.bam");
   system("mv ${input_sampleID}_h1_su.bam ${input_sampleID}_su.bam");
   if($mode eq "standard" || $sequencing_type eq "single-end"){
     system("gunzip -c ${input_sampleID}_soft.fastq.gz >${input_sampleID}_1soft.fastq");
     system("gunzip -c ${input_sampleID}_h1_soft.fastq.gz >>${input_sampleID}_1soft.fastq");
   }
   if($sequencing_type eq "paired-end"){
     system("mv ${input_sampleID}_h1_h1_1.1fq ${input_sampleID}_1.1fq");
     system("mv ${input_sampleID}_h1_h1_2.1fq ${input_sampleID}_2.1fq");
   } else {
     system("mv ${input_sampleID}_h1_h1.1fq ${input_sampleID}.1fq");
   }
 } else {
   $order=1;
   align_to_hg(${input_sampleID},${file_suffix});
   convert_bamtofastq(${input_sampleID});
   if ($data_type eq "RNA-seq") {
     $order=2;
     align_to_hg(${input_sampleID}."_h1",".1fq");
     convert_bamtofastq(${input_sampleID}."_h1");
     system("mv ${input_sampleID}_h1_sm.bam ${input_sampleID}_sm.bam");
     system("mv ${input_sampleID}_h1_su.bam ${input_sampleID}_su.bam");
     if($mode eq "standard" || $sequencing_type eq "single-end"){
       system("gunzip -c ${input_sampleID}_soft.fastq.gz >${input_sampleID}_1soft.fastq");
       system("gunzip -c ${input_sampleID}_h1_soft.fastq.gz >>${input_sampleID}_1soft.fastq");
     } 
     if($sequencing_type eq "paired-end"){
       system("mv ${input_sampleID}_h1_h1_1.1fq ${input_sampleID}_1.1fq");
       system("mv ${input_sampleID}_h1_h1_2.1fq ${input_sampleID}_2.1fq");
     } else {
       system("mv ${input_sampleID}_h1_h1.1fq ${input_sampleID}.1fq");
     }
   } else { 
     if($mode eq "standard" || $sequencing_type eq "single-end"){
       system("gunzip -c ${input_sampleID}_soft.fastq.gz >${input_sampleID}_1soft.fastq");
     }
     if($sequencing_type eq "paired-end"){
       system("mv ${input_sampleID}_h1_1.1fq ${input_sampleID}_1.1fq");
       system("mv ${input_sampleID}_h1_2.1fq ${input_sampleID}_2.1fq");
     } else {
       system("mv ${input_sampleID}_h1.1fq ${input_sampleID}.1fq");
     }
   }
 }

 ##### Filter split reads
 open SF1, "${input_sampleID}_1soft.fastq";
 open SF2, ">${input_sampleID}_1soft.fastq2";
 while (my $tmp1=<SF1>) {
   if ($tmp1=~"^\@soft") {
     my @tmp1=split /\|/, $tmp1;
     if ($tmp1[2]%2 >=1) {
       print SF2 "$tmp1";
       $tmp1=<SF1>;
       print SF2 "$tmp1";
       $tmp1=<SF1>;
       print SF2 "$tmp1";
       $tmp1=<SF1>;
       print SF2 "$tmp1";
     } else {
       $tmp1=<SF1>;
       $tmp1=<SF1>;
       $tmp1=<SF1>;
     }
   }
 }
 close SF1;
 close SF2;
 system ("mv ${input_sampleID}_1soft.fastq2 ${input_sampleID}_1soft.fastq");

 print "~~~~~ step 2\n";
 ##### quality control
 if (defined($QS_cutoff) && $sequencing_type eq "paired-end" && $mode eq "standard"){
   fastq_qc(${input_sampleID},1);
 } elsif (defined($QS_cutoff) && $sequencing_type eq "paired-end" && $mode eq "fast"){
   fastq_qc(${input_sampleID},2);
 } elsif (defined($QS_cutoff) && $sequencing_type eq "single-end"){
   fastq_qc(${input_sampleID},3);
 }
 if(!defined($QS_cutoff)){                        ##### single-end or standard mode without quality control
   if($sequencing_type eq "single-end" || $mode eq "standard"){
     system ("perl ${directory}Scripts/Soft_clipping_filter.pl -length 20 -file ${input_sampleID}_1soft.fastq -o ${input_sampleID}");
     system ("rm ${input_sampleID}_1soft.fastq");
   }
 }

 print "~~~~~ step 3\n";
 ##### align to viral genomes
 if($sequencing_type eq "paired-end"){
   system ("${bwa_d}bwa mem -t $threads -k 19 -r 1.5 -c 100000 -m 50 -T 20 -h 10000 -a -Y -M $virus_genome ${input_sampleID}_1.1fq ${input_sampleID}_2.1fq >${input_sampleID}_vsu.sam");
 } else{	
   system ("${bwa_d}bwa mem -t $threads -k 19 -r 1.5 -c 100000 -m 50 -T 20 -h 10000 -a -Y -M $virus_genome ${input_sampleID}.1fq >${input_sampleID}_vsu.sam");
 }
 system ("${samtools_d}samtools view -bS -@ $thread_1 ${input_sampleID}_vsu.sam >${input_sampleID}_vsu.bam");
 system ("${samtools_d}samtools sort -@ $thread_1 ${input_sampleID}_vsu.bam -o ${input_sampleID}_vsu.sort.bam");
 system ("${samtools_d}samtools view -@ $thread_1 ${input_sampleID}_vsu.sort.bam >${input_sampleID}_vsu.sort.sam");
 system ("rm ${input_sampleID}_vsu.bam");
 system ("rm ${input_sampleID}_vsu.sam");
 system ("touch ${input_sampleID}_all_breakpoint");

 ##### split reads
 if($sequencing_type eq "single-end" || ($sequencing_type eq "paired-end" && $mode eq "standard")){
   system ("${bwa_d}bwa mem -t $threads -k 19 -r 1.5 -c 100000 -m 50 -T 20 -h 10000 -a -Y -M $virus_genome ${input_sampleID}_1sf.fastq >${input_sampleID}_vsoft.sam");
   system ("${samtools_d}samtools view -bS -@ $thread_1 ${input_sampleID}_vsoft.sam >${input_sampleID}_vsoft.bam");
   system ("${samtools_d}samtools sort -@ $thread_1 ${input_sampleID}_vsoft.bam -o ${input_sampleID}_vsoft_sort.bam");
   system ("${samtools_d}samtools view -@ $thread_1 ${input_sampleID}_vsoft_sort.bam >${input_sampleID}_vsoft_sort.sam");
   system ("rm ${input_sampleID}_vsoft.bam");
   system ("perl ${directory}Scripts/Soft_clipping_transfer.pl -f ${input_sampleID}_vsoft_sort.sam -taxonomy $virus_taxonomy -o ${input_sampleID}_vsoft_breakpoint");
   system ("cat ${input_sampleID}_vsoft_breakpoint >>${input_sampleID}_all_breakpoint");
   system ("rm ${input_sampleID}_vsoft_breakpoint");
   system ("rm ${input_sampleID}_vsoft_sort.sam");
   system ("rm ${input_sampleID}_vsoft.sam");
 }

 print "~~~~~ step 4\n";
 if($sequencing_type eq "paired-end"){ 
   create_type();
   system ("perl ${directory}Scripts/Break_point_calling.pl -taxonomy $virus_taxonomy -type ${input_sampleID}.type -position ${input_sampleID}_sm.sam -virus ${input_sampleID}_vsu.sort.sam -o ${input_sampleID}");
   system ("rm ${input_sampleID}_sm.sam");
   system ("cat ${input_sampleID}_breakpoint >>${input_sampleID}_all_breakpoint");
   system ("rm ${input_sampleID}_breakpoint");
 }
 system ("perl ${directory}Scripts/Reads_summary.pl -i ${input_sampleID}_all_breakpoint -o ${input_sampleID}_summary");
 system ("rm ${input_sampleID}_all_breakpoint");
 system ("perl ${directory}Scripts/Order_by_virus_sequence.pl -i ${input_sampleID}_summary -o ${input_sampleID}.integration");
 my $double_length_insertsize=$length_insertsize * 2;
 system ("perl ${directory}Scripts/Filtered_single_reads.pl -i ${input_sampleID}.integration -o ${input_sampleID}.integration2 -r $double_length_insertsize");
 system ("mv ${input_sampleID}.integration2 ${input_sampleID}.integration");
 print "~~~~~ step 5\n";
 ##### Infection screening
 if($mode eq "standard"){
   detect_virus();
   system ("rm ${input_sampleID}_vsu.sort.sam");
 }

 print "~~~~~ step 6\n";
 ##### obtain fuq reads;
 system ("perl ${directory}Scripts/Filtered_fastq.pl $input_sampleID");
 if ($sequencing_type eq "single-end"){
   system ("rm ${input_sampleID}_1.1fuq ${input_sampleID}_2.1fuq");
 }

 ##### reciprocal alignment
 if(defined($align_back_to_human)){
   reciprocal_alignment();
 }

 print "~~~~~ step 7\n";
 ##### repeats
 if(defined($repeat)){
   filter_repeats();
 }
 system ("perl ${directory}Scripts/Filtered_integration-HBV_final.pl ${input_sampleID} >${input_sampleID}.21");
 system ("perl ${directory}Scripts/Remove_redundancy_split_reads.pl ${input_sampleID}");
 system ("rm ${input_sampleID}.21");
 system ("mv ${input_sampleID}.2 ${input_sampleID}.3");

 print "~~~~~ step 8\n";
 #####  output
 system ("perl ${directory}Scripts/Assign_reads-filter.pl -i ${input_sampleID}.3 -o ${input_sampleID}.f -r $double_length_insertsize -m f");
 system ("perl ${directory}Scripts/Result_filtered.pl ${input_sampleID}.f >${input_sampleID}_f");
 system ("perl ${directory}Scripts/Result_visual3-3.pl ${input_sampleID}_f ${input_sampleID} 12 80 >${input_sampleID}_f2");
 system ("perl ${directory}Scripts/Result_finalize3.pl ${input_sampleID}_f2 >${input_sampleID}_f22");
 system ("perl ${directory}Scripts/Results_get.pl ${input_sampleID}_f22 >${input_sampleID}.virus_f");
 my $cmd="";
 if ($sequencing_type eq "single-end"){
   $cmd = q(awk '{if(($15>=50)&&(($39>2 && $41>=30) || ($39==2 && $41>=50)))print$0}');
 } else {
   $cmd = q(awk '{if(($15>=50)&&($24>0)&&(($39>2 && $41>=30) || ($39==2 && $41>=50)))print$0}');
 }
 system ("$cmd ${input_sampleID}.virus_f >${input_sampleID}.virus_f2");
 system ("perl ${directory}Scripts/Extract_specific_loci_final_visualization.pl ${input_sampleID}_f2 ${input_sampleID}.virus_f2 >${input_sampleID}.visualization");
 system ("perl ${directory}Scripts/Fine_mapped_VI.pl ${input_sampleID}.visualization >${input_sampleID}.fine_mapped");

 ##### align to vector db
 if($sequencing_type eq "paired-end"){
   system ("${bwa_d}bwa mem -t $threads -k 19 -r 1.5 -c 100000 -m 50 -T 20 -h 10000 -a -Y -M $vector_db ${input_sampleID}_1.1fuq ${input_sampleID}_2.1fuq >${input_sampleID}_vector.sam");
   system ("perl ${directory}Scripts/Filter_vector.pl ${input_sampleID} vector");
 }
 if ($sequencing_type eq "single-end" || $mode eq "standard") {
   system ("${bwa_d}bwa mem -t $threads -k 19 -r 1.5 -c 100000 -m 50 -T 20 -h 10000 -a -Y -M $vector_db ${input_sampleID}_1sf.fuq >${input_sampleID}_vector_sf.sam");
   system ("perl ${directory}Scripts/Filter_vector.pl ${input_sampleID} vector_sf");
 } 
 my @path_current2=split /\//, ${input_sampleID};
 my $path_current2;
 if($#path_current2>0){
   my $path_current2=join("/",@path_current2[0..$#path_current2-1]);
   my $input3=$path_current2[$#path_current2];
 }
 else {
   $path_current2= Cwd::cwd();
   $path_current2=$path_current2."/";
 }
 ##### delete files
 system ("rm ${input_sampleID}.integration");
 system ("rm ${input_sampleID}.fasta");
 system ("rm ${input_sampleID}_summary");
 system ("rm ${input_sampleID}.f");
 system ("rm ${input_sampleID}_f");
 system ("rm ${input_sampleID}_f22");
 if(-e ${input_sampleID}."_1sf.others"){
   system ("rm ${input_sampleID}_1sf.others");
 }
 if($sequencing_type eq "single-end" && $mode eq "fast"){
   system ("rm ${input_sampleID}.1fq");
   system ("rm ${input_sampleID}_vsu.sort.bam");
   system ("rm ${input_sampleID}_vsu.sort.sam");
 }

##### validation function
} elsif ($function eq "validate"){
  if (!(-d ${directory}."Database/GI")) {
    mkdir ${directory}."Database/GI/";
  }
  @groups=("target","vicaller","human");
  v_obtain_seq();
  v_index_GI();
  v_bwa_validate();
  v_combine();
  v_blat_validate();
  v_blastn_validate();
  v_convert();
  v_summary();

##### calculate function
} elsif ($function eq "calculate"){
  my $position1=$Position-10000;
  my $position2=$Position+10000;

  if (${File_suffix_bam} =~ "fq" || ${File_suffix_bam} =~ "fastq") {
    system ("${samtools_d}samtools sort -@ $thread_1 ${input_sampleID}_h.bam -o ${input_sampleID}_s_h.bam");
    system ("${samtools_d}samtools index ${input_sampleID}_s_h.bam");
    system ("${samtools_d}samtools view ${input_sampleID}_s_h.bam ${Chr}:${position1}-${position2} >${input_sampleID}_${Chr}_${Position}_h.sort.sam");
  } elsif (${File_suffix_bam} =~ ".bam" && defined($Index_sort)){
    system ("${samtools_d}samtools view ${input_sampleID}${File_suffix_bam} ${Chr}:${position1}-${position2} >${input_sampleID}_${Chr}_${Position}_h.sort.sam");
  } elsif (${File_suffix_bam} =~ ".bam" && !defined($Index_sort)) {
    system ("${samtools_d}samtools sort -@ $thread_1 ${input_sampleID}${File_suffix_bam} -o ${input_sampleID}_s${File_suffix_bam}");
    system ("${samtools_d}samtools index ${input_sampleID}_s${File_suffix_bam}");
    system ("${samtools_d}samtools view ${input_sampleID}_s${File_suffix_bam} ${Chr}:${position1}-${position2} >${input_sampleID}_${Chr}_${Position}_h.sort.sam");
         }
  system ("sort ${input_sampleID}_${Chr}_${Position}_h.sort.sam | uniq > ${input_sampleID}_${Chr}_${Position}_h.sort.sam2"); 
  system ("perl ${directory}Scripts/C_filtered1.pl ${input_sampleID}_${Chr}_${Position}_h.sort.sam2 >${input_sampleID}_${Chr}_${Position}_bp1");
  system ("perl ${directory}Scripts/C_filtered2.pl ${input_sampleID}_${Chr}_${Position}_bp1 >${input_sampleID}_${Chr}_${Position}_bp2");
  system ("perl ${directory}Scripts/C_filtered3.pl ${input_sampleID}_${Chr}_${Position}_bp2 $Chr $Position >${input_sampleID}_${Chr}_${Position}.bp");
  system ("rm ${input_sampleID}_${Chr}_${Position}_bp2 ${input_sampleID}_${Chr}_${Position}_bp1");
  system ("rm ${input_sampleID}_${Chr}_${Position}_h.sort.sam ${input_sampleID}_${Chr}_${Position}_h.sort.sam2");
  system ("perl ${directory}Scripts/C_allele_fraction.pl ${input_sampleID}_${Chr}_${Position} ${input_sampleID} $Chr $Position $Number_reads $Breakpoint");
}

##### output
##### fine_mapped file
 my %output=();
 my $output_header="";
 my @output_line=();
 my $header="Sample_ID  VIcaller_mode   QC      Reciprocal_alignment    Candidate_virus GI      Chr.    Start   End     No._chimeric_reads      No._split_reads Upstream_breakpoint_on_human    Downstream_breakpoint_on_human  Upstream_breakpoint_on_virus    Downstream_breakpoint_on_virus  Information_of_both_upstream_and_downstream_breakpoints Integration_site_in_the_human_genome    Integration_allele_fraction     No._reads_supporting_nonVI      No._reads_supporting_VI Average_alignment_score Is_cell_line_contamination Is_vector Validation_chimeric_confident   Validation_chimeric_weak        Validation_chimeric_false       Validation_split_confident      Validation_split_weak   Validation_split_false\n";
if ($function eq "detect"){
 print "~~~~~ step 9 output of VIcaller detect function\n";

 ##### open file_mapped
 open FINE, "${input_sampleID}.fine_mapped";
 while (<FINE>){
   @output_line=split;
   $output_header=join ("_",(@output_line[0..3],$output_line[5]));
   @{$output{$output_header}} = ("-") x 29;
   $output{$output_header}[0]=$output_line[0];
   if($mode eq "standard"){
     $output{$output_header}[1] = "standard";
   } elsif ($mode eq "fast"){
     $output{$output_header}[1] = "fast";
   }
   if(defined($QS_cutoff)){
     $output{$output_header}[2] = "yes";
   } else {
     $output{$output_header}[2] = "no";
   }
   if(defined($align_back_to_human)){
     $output{$output_header}[3] = "yes";
   } else {
     $output{$output_header}[3] = "no";
   }
   @{$output{$output_header}}[4..5]=@output_line[4..5];
   @{$output{$output_header}}[6..8]=@output_line[1..3];
   if(($output_line[8] ne "-" && $output_line[17] eq "-") || ($output_line[8] ne "-" && $output_line[17] ne "-" && $output_line[16]<=$output_line[7])){
     if ($output_line[7] eq "-") {
       $output{$output_header}[15]="E(++);";
       $output{$output_header}[11]=$output_line[9];
       $output{$output_header}[14]=$output_line[12];
     } else {
       $output{$output_header}[15]="D(++);";
       $output{$output_header}[11]=$output_line[7];
       $output{$output_header}[14]=$output_line[10];
     }   
   } elsif (($output_line[8] eq "-" && $output_line[17] ne "-") || ($output_line[8] ne "-" && $output_line[17] ne "-" && $output_line[16]>=$output_line[7])){
     if ($output_line[16] eq "-") {
       $output{$output_header}[15]="E(+-);";
       $output{$output_header}[11]=$output_line[18];
       $output{$output_header}[13]=$output_line[20];
     } else {
       $output{$output_header}[15]="D(+-);";
       $output{$output_header}[11]=$output_line[16];
       $output{$output_header}[13]=$output_line[19];
     }
   } elsif ($output_line[8] eq "-" && $output_line[17] eq "-") {
     $output{$output_header}[15]="na;";
   }
  
   if(($output_line[26] ne "-" && $output_line[35] eq "-") || ($output_line[26] ne "-" && $output_line[35] ne "-" && $output_line[34]<=$output_line[25])){
     if ($output_line[25] eq "-") {
       $output{$output_header}[15]=$output{$output_header}[15]."e(-+)";
       $output{$output_header}[12]=$output_line[26];
       $output{$output_header}[14]=$output_line[30];
     } else {
       $output{$output_header}[15]=$output{$output_header}[15]."d(-+)";
       $output{$output_header}[12]=$output_line[25];
       $output{$output_header}[14]=$output_line[28];
     }
   } elsif (($output_line[26] eq "-" && $output_line[35] ne "-") || ($output_line[26] ne "-" && $output_line[35] ne "-" && $output_line[34]>=$output_line[25])){
     if ($output_line[34] eq "-") {
       $output{$output_header}[15]=$output{$output_header}[15]."e(--)";
       $output{$output_header}[12]=$output_line[35];
       $output{$output_header}[13]=$output_line[38];
     } else {
       $output{$output_header}[15]=$output{$output_header}[15]."d(--)";
       $output{$output_header}[12]=$output_line[34];
       $output{$output_header}[13]=$output_line[37];
     }
   } elsif ($output_line[26] eq "-" && $output_line[35] eq "-") {
     $output{$output_header}[15]=$output{$output_header}[15]."na";
   }

   if ($output{$output_header}[11] ne "-" && $output{$output_header}[12] ne "-"){
     if (($output{$output_header}[15] =~ "D" && $output{$output_header}[15] =~ "d") || ($output{$output_header}[15] =~ "E" && $output{$output_header}[15] =~ "e")){
       $output{$output_header}[16]=int(($output{$output_header}[12]+$output{$output_header}[11])/2);
     } elsif ($output{$output_header}[15] =~ "D"){
       $output{$output_header}[16]=$output{$output_header}[11];
     } elsif ($output{$output_header}[15] =~ "d") {
       $output{$output_header}[16]=$output{$output_header}[12];
     }
   } elsif ($output{$output_header}[11] ne "-" && $output{$output_header}[12] eq "-") {
     $output{$output_header}[16]=$output{$output_header}[11];
   } elsif ($output{$output_header}[12] ne "-" && $output{$output_header}[11] eq "-") {
     $output{$output_header}[16]=$output{$output_header}[12];
   }
   $output{$output_header}[15] =~ s/e/E/;
   $output{$output_header}[15] =~ s/d/D/;

   ##### read cell line
   my %cell_line=();
   my @cell_line=();
   my $col_ID=0;
   if ($build eq "hg19"){
     $col_ID=3;
   } else {
     $col_ID=4;
   }
   open CELL,"$cell_line";
   while(<CELL>){
     @cell_line=split;
     my $cell1=$cell_line[$col_ID];
     $cell1=~s/\:/ /;
     $cell1=~s/\-/ /;
     my @cell1=split /\s+/, $cell1;
     if ($output{$output_header}[6] =~ "Chr"){
       $cell1[0] =~ s/chr/Chr/;
     } elsif ($output{$output_header}[6] =~ "chr") {
       $cell1[0] =~ s/chr/chr/;
     } elsif (!($output{$output_header}[6] =~ "chr") && !($output{$output_header}[6] =~ "Chr")) {
       $cell1[0] =~ s/chr//;
     }
     @{$cell_line{$cell_line[$col_ID]}}=(@cell_line[0..2],@cell1);
   }
   close CELL;
   my @key_cell1 = keys %cell_line;
   my $Is_cell="-";
   for (my $ci=0;$ci<@key_cell1;$ci++){
     my @cell2 = @{$cell_line{$key_cell1[$ci]}};
     if ($cell2[3] eq $output{$output_header}[6]){
       if (abs($cell2[4] - $output{$output_header}[7])<=1000 || abs($cell2[4] - $output{$output_header}[8])<=1000 || ($cell2[4]>=$output{$output_header}[7] && $cell2[4]<=$output{$output_header}[8])){ $output{$output_header}[21]="Yes";}
       if ($cell2[5] ne "N" && $cell2[5] ne "NA"){
         if(abs($cell2[5] - $output{$output_header}[7])<=1000 || abs($cell2[5] - $output{$output_header}[8])<=1000 || ($cell2[5]>=$output{$output_header}[7] && $cell2[5]<=$output{$output_header}[8])){$output{$output_header}[21]="Yes";}
       }
     }
   }
 }

 ##### open vector file
 my %vec_line=();
 my @vec_line=();
 my %v_chimeric=();
 my %v_split=();
 system ("cat ${input_sampleID}_vector.hmap ${input_sampleID}_vector_sf.hmap > ${input_sampleID}_vector.hmapped");
 open VEC, "${input_sampleID}_vector.hmapped";
 while (<VEC>){
   @vec_line=split;
   if (exists ($vec_line{$vec_line[1]})){
     if ($vec_line{$vec_line[1]} eq "P") {
       next;
     } elsif (($vec_line{$vec_line[1]} eq "L" && $vec_line[0] eq "L") || ($vec_line{$vec_line[1]} eq "R" && $vec_line[0] eq "R")) {
       next;
     } elsif ($vec_line{$vec_line[1]} eq "L" && $vec_line[0] eq "R" && (($vec_line[7]-$vec_line[6]+1) >= 50 || (($vec_line[7]-$vec_line[6]+1)>=30 && ($vec_line[7]-$vec_line[6]+1)/$vec_line[8]>=0.5))) {
       $vec_line{$vec_line[1]}="P";
     } elsif ($vec_line{$vec_line[1]} eq "R" && $vec_line[0] eq "L" && (($vec_line[4]-$vec_line[3]+1) >= 50 || (($vec_line[4]-$vec_line[3]+1)>=30 && ($vec_line[4]-$vec_line[3]+1)/$vec_line[5]>=0.5))) {
       $vec_line{$vec_line[1]}="P";
     }
   } else {
     if (($vec_line[7]-$vec_line[6]+1) >= 50 || ($vec_line[7]-$vec_line[6]+1) eq $vec_line[8] || (($vec_line[7]-$vec_line[6]+1)>=30 && ($vec_line[7]-$vec_line[6]+1)/$vec_line[8]>=0.5)) {
       $vec_line{$vec_line[1]}="R";
     }
     if (($vec_line[4]-$vec_line[3]+1) >= 50 || ($vec_line[4]-$vec_line[3]+1) eq $vec_line[5] || (($vec_line[4]-$vec_line[3]+1)>=30 && ($vec_line[4]-$vec_line[3]+1)/$vec_line[5]>=0.5)) {
       if (exists ($vec_line{$vec_line[1]})){
         $vec_line{$vec_line[1]}="P";
       } else {
         $vec_line{$vec_line[1]}="L";
       }
     }
   }
 }

 ##### read virus_f2 file
 my @f2_line=();
 my $f2_header="";
 open F2,"${input_sampleID}.virus_f2";
 while (my $f2_line=<F2>){
   @f2_line=split /\s+/, $f2_line;
   ##### check vector
   $f2_line[0]=~s/ad://;
   my $f2_header2 = $f2_line[0]."_".$f2_line[15]."_".$f2_line[35]."_".$f2_line[30]."_".$f2_line[31]."_".$f2_line[32]."_".$f2_line[33];
   my %f2_header2=();
   $f2_header2{$f2_header2}=0;
   open VIS,"${input_sampleID}.visualization";
   while (my $v_line_temp1=<VIS>){
     my @f2_title_O1=();
     my $f2_title3="";
     my $v_ok=0;
     my @v_line=split /\s+/, $v_line_temp1;
     #####  print "VIS: @v_line\n";
     unless (@v_line) {next;}
     if($v_line[0] eq "O1"){
       $v_ok=0;
       my $v_title2=$v_line[2]."_".$v_line[1]."_".$v_line[22]."_".$v_line[17]."_".$v_line[18]."_".$v_line[19]."_".$v_line[20];
       $f2_title3="";
       $v_line[21]= lc $v_line[21];
       if($v_line[21] =~ "unknown") {
         $v_line[21]=$v_line[22];
       }
       if (exists ($f2_header2{$v_title2})) {
         $v_ok=1;
         @f2_title_O1=();
         @f2_title_O1=@v_line;
         $f2_title3=$v_line[1]." ".$v_line[18]." ".$v_line[19]." ".$v_line[20]." ".$f2_header2{$v_title2}." ".$v_line[22];
       }
     } elsif ($v_line[0] eq "O2" && $v_ok==1) {
       if (exists ($vec_line{$v_line[1]})){
         if ($v_line[4] eq "S1"){
           if (exists ($v_split{$v_line[0]})){
             $v_split{$v_line[0]}="Vector";
           } else {
             $v_split{$v_line[0]}="Vector";
           }
         } elsif ($v_line[4] eq "S") {
           if ($vec_line{$v_line[1]} eq "P" || $vec_line{$v_line[1]} eq "R") {
             $v_chimeric{$v_line[0]}="Vector";
           } else {
             $v_chimeric{$v_line[0]}="Non-vector";
           }
         } elsif ($v_line[4] eq "F") {
           if ($vec_line{$v_line[1]} eq "P" || $vec_line{$v_line[1]} eq "L") {
             $v_chimeric{$v_line[0]}="Vector";
           } else {
             $v_chimeric{$v_line[0]}="Non-vector";
           }
         }
       } else {
         if ($v_line[4] eq "S1") {
           $v_split{$v_line[0]}="Non-vector";
         } else {
           $v_chimeric{$v_line[0]}="Non-vector";
         }
       }
     }
   }
   close VIS;

   ##### vector organize
   my @key_chimeric=();
   my @key_split=();
   my $v_chimeric_count=0;
   my $v_split_count=0;
   for(my $kc=0;$kc<@key_chimeric;$kc++){
     if ($v_chimeric{$key_chimeric[$kc]} eq "Vector"){$v_chimeric_count++;}
   }
   for (my $ks=0;$ks<@key_split;$ks++){
     if ($v_split{$key_split[$ks]} eq "Vector"){$v_split_count++;}
   }

   ##### organize;
   $f2_header=join ("_",($f2_line[15],@f2_line[31..33],$f2_line[35]));
   if(exists($output{$f2_header})){
     $output{$f2_header}[9]=$#key_chimeric+1;
     $output{$f2_header}[10]=$#key_split+1;
     $output{$f2_header}[19]=$f2_line[38];
     $output{$f2_header}[20]=$f2_line[40];
     if (%v_chimeric && %v_split) {
       if ($v_chimeric_count eq ($#key_chimeric+1) && $v_split_count eq ($#key_split+1)) {
         $output{$f2_header}[22] = "Yes";
       }
     } elsif (%v_chimeric) {
       if ($v_chimeric_count eq ($#key_chimeric+1)){
         $output{$f2_header}[22] = "Yes";}
     } elsif (%v_split) {
       if ($v_split_count eq ($#key_split+1)) {
         $output{$f2_header}[22] = "Yes";
       }
     } else {
       $output{$f2_header}[22] = "-";
     }
   }
 }

  ##### output 
  my @list=sort{$output{$b}->[19]<=>$output{$a}->[19]} keys %output;
  open OUT, ">${input_sampleID}.output";
  print OUT "$header";
  $" = "\t";
  for(my $i=0;$i<@list;$i++){
    print OUT "@{$output{$list[$i]}}\n";
  }
  close OUT; 
}

##### read validation file
 my $cs3_header="";
 my @cs3_line=();
 %output=();
 my @out_line="";
 my $out_header="";
 if ($function eq "validate" && -e ${input_sampleID}.".output") {
   print "~~~~~ step 9 output of VIcaller validate function\n";
   open OUT, "${input_sampleID}.output";
   my $out_temp=<OUT>;
   while (<OUT>){
     @out_line=split;
     $out_header=join ("_",($out_line[0],@out_line[6..8],$out_line[5]));
     @{$output{$out_header}}=@out_line;
   }

   open CS3, "${String}.CS3";
   while (<CS3>){
     @cs3_line=split;
     $cs3_header=join ("_",($cs3_line[57],$cs3_line[58],$cs3_line[59],$cs3_line[60],$cs3_line[62]));
     if(exists($output{$cs3_header})){
       @{$output{$cs3_header}}[23..25]=@cs3_line[2..4];
       @{$output{$cs3_header}}[26..28]=@cs3_line[6..8];
     }
   }
   my @list=sort{$output{$b}->[19]<=>$output{$a}->[19]} keys %output;
   open OUT, ">${input_sampleID}.output";
   print OUT "$header";
   $" = "\t";
   for(my $i=0;$i<@list;$i++){
     print OUT "@{$output{$list[$i]}}\n";
   }
   close OUT;
 }

##### read bp2 file
 if ($function eq "calculate" && -e ${input_sampleID}."_".${Chr}."_".${Position}.".bp2"){
   print "~~~~~ step 9 output of VIcaller calculate function\n";
   my @out_line="";
   my $out_header="";
   my @bp2_line=();
   %output=();
   open BP2,"${input_sampleID}_${Chr}_${Position}.bp2";
   my $bp2_line=<BP2>;
   $bp2_line=<BP2>;
   @bp2_line= split /\s+/, $bp2_line;
   if (-e ${input_sampleID}.".output") {
     open OUT, "${input_sampleID}.output";
     my $out_temp=<OUT>;
     while (<OUT>){
       @out_line=split;
       $out_header=$out_line[0]."_".$out_line[6]."_".$out_line[16];
       @{$output{$out_header}}=@out_line;
     }
     close OUT;
     open OUT,">${input_sampleID}.output";
     print OUT "$header";
     $" = "\t";
     my @list=sort{$output{$b}->[19]<=>$output{$a}->[19]} keys %output;
     for(my $i=0;$i<@list;$i++){
       if ($list[$i] eq $bp2_line[0]){
         ${$output{$list[$i]}}[17]=$bp2_line[3];
         ${$output{$list[$i]}}[18]=$bp2_line[1];
       }
       print OUT "@{$output{$list[$i]}}\n";
     }
     close OUT;
   }
 }

##### get local time
@months = qw( Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec );
@days = qw(Sun Mon Tue Wed Thu Fri Sat Sun);
($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime();
print "\n\n\n~~~~~ Done\n#############################################################\n";
print "###############                              \n";
print "###############      $hour:$min:$sec, $months[$mon].$mday\n";
print "###############      Author: Xun Chen (Xun.Chen\@uvm.edu and xunchen85\@gmail.com)\n";
print "###############                              \n";
print "#############################################################\n\n\n\n";

##### all sub functions
##### detect function
##### detect virus
sub detect_virus {
  system("perl ${directory}Scripts/Infection_filtered_detection-1-092717.pl ${input_sampleID}_vsu.sort.sam ${input_sampleID}.type 30 >${input_sampleID}_In_infection2");
  system("perl ${directory}Scripts/Infection_filtered_detection-2-092717.pl ${input_sampleID}_In_infection2 $virus_list >${input_sampleID}_In_infection22");
  system("perl ${directory}Scripts/Infection_filtered_detection-3-092717.pl ${input_sampleID}_In_infection22 >${input_sampleID}_In.virus_v");
  system("perl ${directory}Scripts/Infection_filtered_detection-4-092717.pl ${input_sampleID}_In_infection22 >${input_sampleID}_In.virus_g");
  system("perl ${directory}Scripts/Infection_Summary_virus-1-092717.pl ${input_sampleID}_In.virus_v ${input_sampleID}_In.virus_g $virus_taxonomy $virus_list >${input_sampleID}_In.virus_v2");
  system("perl ${directory}Scripts/Infection_Summary_virus-2-092717.pl ${input_sampleID}_In.virus_v2 5 >${input_sampleID}_In.virus_v3");
  if (-e ${input_sampleID}."_In_infection2") {system("rm ${input_sampleID}_In_infection2");}
  if (-e ${input_sampleID}."_In_infection22") {system("rm ${input_sampleID}_In_infection22");}
  if (-e ${input_sampleID}."_In.virus_v") {system("rm ${input_sampleID}_In.virus_v");}
  if (-e ${input_sampleID}."_In.virus_g") {system("rm ${input_sampleID}_In.virus_g");}
  if (-e ${input_sampleID}."_In.virus_v2") {system("rm ${input_sampleID}_In.virus_v2");}
}

##### reciprocal alignment
sub reciprocal_alignment {
  if($sequencing_type eq "paired-end"){
    system ("${bwa_d}bwa mem -t $threads -k 19 -r 1.5 -c 100000 -m 50 -T 20 -h 10000 -a -Y -M $human_genome ${input_sampleID}_1.1fuq ${input_sampleID}_2.1fuq >${input_sampleID}_hpe.sam");
    system ("perl ${directory}Scripts/Realign_filter2.pl $input_sampleID");
    system ("sort -k 2 -k 3 -k 4 -k 5 -k 6 -k 7 -k 8 -k 9 -k 10 ${input_sampleID}.hmap |uniq >${input_sampleID}.hmapped");
    if (-e ${input_sampleID}.".hmap") {system ("rm ${input_sampleID}.hmap");}
    if (-e ${input_sampleID}.".hpe.sam") {system ("rm ${input_sampleID}_hpe.sam");}
  }
}

##### check repeat
sub filter_repeats {
  my $path_current= Cwd::cwd();
  my @directory1=split /\//, ${input_sampleID};
  if($#directory1>0){
    my $directory1=join("/",@directory1[0..$#directory1-1]);
    my $input2=$directory1[$#directory1];
    print "$path_current\n";
    chdir $path_current."/".$directory1;
    my $path_current1= Cwd::cwd();
    print "$path_current1\n";
    system ("${repeatmasker_d}trf ${input2}.fasta 2 5 7 80 10 50 2000 -h -f -d -m");
    system ("${dust_d}bin/dust ${input2}.fasta.2.5.7.80.10.50.2000.mask >${input2}.fasta2");
    if (-e ${input2}.".fasta.2.5.7.80.10.50.2000.mask") {system("rm ${input2}.fasta.2.5.7.80.10.50.2000.mask");}
    system ("perl ${directory}Scripts/Exclude_trf_dust_reads061617.pl ${input2}.fasta2 >${input2}.fasta3");
    system ("${repeatmasker_d}RepeatMasker -e rmblast -pa $threads -species human -s ${input2}.fasta3");
    system ("perl ${directory}Scripts/Reformat_repeatmasker_result.pl ${input2} >${input2}.repeat");
    system ("perl ${directory}Scripts/Reformat_trf_dust_result2.pl ${input2} >${input2}.repeat2");
    chdir $path_current;
  } else {
    system ("${repeatmasker_d}trf ${input_sampleID}.fasta 2 5 7 80 10 50 2000 -h -f -d -m");
    system ("${dust_d}bin/dust ${input_sampleID}.fasta.2.5.7.80.10.50.2000.mask >${input_sampleID}.fasta2");
    if (-e ${input_sampleID}.".fasta.2.5.7.80.10.50.2000.mask") {system ("rm ${input_sampleID}.fasta.2.5.7.80.10.50.2000.mask");}
    system ("perl ${directory}Scripts/Exclude_trf_dust_reads061617.pl ${input_sampleID}.fasta2 >${input_sampleID}.fasta3");
    system ("${repeatmasker_d}RepeatMasker -e rmblast -pa $threads -species human -s ${input_sampleID}.fasta3");
    system ("perl ${directory}Scripts/Reformat_repeatmasker_result.pl ${input_sampleID} >${input_sampleID}.repeat");
    system ("perl ${directory}Scripts/Reformat_trf_dust_result2.pl ${input_sampleID} >${input_sampleID}.repeat2");
  }
  if (-e ${input_sampleID}.".repeat") {system ("rm ${input_sampleID}.repeat");}
  if (-e ${input_sampleID}.".fasta3") {system ("rm ${input_sampleID}.fasta3");}
  if (-e ${input_sampleID}.".fasta3.cat") {system ("rm ${input_sampleID}.fasta3.cat");}
  if (-e ${input_sampleID}.".fasta3.masked") {system ("rm ${input_sampleID}.fasta3.masked");}
  if (-e ${input_sampleID}.".fasta3.out") {system ("rm ${input_sampleID}.fasta3.out");}
  if (-e ${input_sampleID}.".fasta3.tbl") {system ("rm ${input_sampleID}.fasta3.tbl");}
  system ("rm -rf ${input_sampleID}*RMoutput");
  if (-e ${input_sampleID}.".fasta.2.5.7.80.10.50.2000.dat") {system ("rm ${input_sampleID}.fasta.2.5.7.80.10.50.2000.dat");}
  if (-e ${input_sampleID}.".fasta2") {system ("rm ${input_sampleID}.fasta2");}
}

##### create type file
sub create_type { 
  system ("${samtools_d}samtools view ${input_sampleID}_sm.bam >${input_sampleID}_sm.sam"); 
  open SM,"${input_sampleID}_sm.sam";
  open TYPE,">${input_sampleID}.type";
  my $sm_1="";
  while($sm_1=<SM>){
    my @sm_1=split /\s+/, $sm_1;
    if($sm_1[1]%256>=128) {
      print TYPE "$sm_1[0] L\n";
    } else {
      print TYPE "$sm_1[0] R\n";
    }
  }
  close TYPE;
  close SM;
}

##### quality control main
sub fastq_qc {
  my ($input1,$seq_type1) = @_;
  if($seq_type1 eq 1){
    fastq_split_qc($input1);
    fastq_PE_qc($input1);
  } elsif ($seq_type1 eq 2){
    fastq_PE_qc($input1);
  } elsif ($seq_type1 eq 3){            
    fastq_split_qc($input1);
    if ($mode eq "standard") {
      system ("perl ${NGSQCToolkit_d}Trimming/TrimmingReads_sanger.pl -i ${input1}.1fq -q $QS_cutoff -n 30");
      system ("perl ${NGSQCToolkit_d}QC/IlluQC_PRLL.pl -se ${input1}.1fq_trimmed 2 1 -c $threads -l 80 -s $QS_cutoff");
      my @directory1=split /\//, ${input1};
      my $directory1=join("/",@directory1[0..$#directory1-1]);
      my $input2=$directory1[$#directory1];
      if($#directory1>0){
        system ("mv ${directory1}/IlluQC_Filtered_files/${input2}.1fq_trimmed_filtered ${input1}.1fq");
      } else {
        system ("mv IlluQC_Filtered_files/${input1}.1fq_trimmed_filtered ${input1}.1fq");
      }
      system ("rm ${input1}.1fq_trimmed");
    }
  }
}

##### fastq quality control for paired_end_reads
sub fastq_PE_qc {
  my ($input1)=@_;
  system ("perl ${NGSQCToolkit_d}Trimming/TrimmingReads_sanger.pl -i ${input1}_1.1fq -irev ${input1}_2.1fq -q $QS_cutoff -n 30");
  system ("perl ${NGSQCToolkit_d}QC/IlluQC_PRLL.pl -pe ${input1}_1.1fq_trimmed ${input1}_2.1fq_trimmed 2 1 -c $threads -l 80 -s $QS_cutoff");
  my @directory1=split /\//, ${input1};
  my $directory1=join("/",@directory1[0..$#directory1-1]);
  my $input2=$directory1[$#directory1];
  if($#directory1>0){
    system ("mv ${directory1}/IlluQC_Filtered_files/${input2}_1.1fq_trimmed_filtered ${input1}_1.1fastq");
    system ("mv ${directory1}/IlluQC_Filtered_files/${input2}_2.1fq_trimmed_filtered ${input1}_2.1fastq");
  } else {
    system ("mv IlluQC_Filtered_files/${input1}_1.1fq_trimmed_filtered ${input1}_1.1fastq");
    system ("mv IlluQC_Filtered_files/${input1}_2.1fq_trimmed_filtered ${input1}_2.1fastq");
  }
  system ("echo ${input1}_1.1fastq >${input1}_list");
  system ("echo ${input1}_2.1fastq >>${input1}_list");
  system ("${fastuniq_d}source/fastuniq -i ${input1}_list -t q -o ${input1}_1.1fq -p ${input1}_2.1fq");
  system ("rm ${input1}_list");
  system ("rm ${input1}_1.1fq_trimmed");
  system ("rm ${input1}_2.1fq_trimmed");
  system ("rm ${input1}_1.1fastq");
  system ("rm ${input1}_2.1fastq");
}

##### fastq quality control for split reads
sub fastq_split_qc {
  my ($input1) = @_;
  system ("perl ${NGSQCToolkit_d}Trimming/TrimmingReads_sanger.pl -i ${input1}_1soft.fastq -q $QS_cutoff");
  system ("perl ${NGSQCToolkit_d}QC/IlluQC_PRLL.pl -se ${input1}_1soft.fastq_trimmed 2 1 -c $threads -l 80 -s $QS_cutoff");
  my @directory1=split /\//, ${input1};
  my $directory1=join("/",@directory1[0..$#directory1-1]);
  my $input2=$directory1[$#directory1];
  if($#directory1>0) {
    system ("mv ${directory1}/IlluQC_Filtered_files/${input2}_1soft.fastq_trimmed_filtered ${input1}_1soft.fastq");
  } else {
    system ("mv IlluQC_Filtered_files/${input1}_1soft.fastq_trimmed_filtered ${input1}_1soft.fastq");
  }
  system ("rm ${input1}_1soft.fastq_trimmed");
  system ("perl ${directory}Scripts/Soft_clipping_filter.pl -length 20 -file ${input_sampleID}_1soft.fastq -o ${input_sampleID}");
  system ("rm ${input_sampleID}_1soft.fastq");
}

##### align to the human reference genome
sub align_to_hg {
  my ($input1,$suffix1)=@_;
  if ($data_type eq "RNA-seq" && $order ne 2){
    my $suffix2=$suffix1;
    $suffix2=~s/.gz$//;
    if($sequencing_type eq "paired-end"){
      if($suffix1 =~ "gz"){
        system ("gunzip -c ${input1}_1${suffix1} >${input1}_1${suffix2}");
        system ("gunzip -c ${input1}_2${suffix1} >${input1}_2${suffix2}");
      }
      system ("${tophat_d}tophat -p $threads -o ${input1} $human_genome_tophat ${input1}_1${suffix2} ${input1}_2${suffix2}");
    } else {
      if($suffix1 =~ "gz"){
        system ("gunzip -c ${input1}${suffix1} >${input1}${suffix2}");
      }
      system ("${tophat_d}tophat -p $threads -o ${input1} $human_genome_tophat ${input1}${suffix2}");
    }
    system ("${samtools_d}samtools merge -f ${input1}_h.bam ${input1}/accepted_hits.bam ${input1}/unmapped.bam");
    system ("${samtools_d}samtools view ${input1}_h.bam >${input1}_h.sam"); 
  } else {
    if($sequencing_type eq "paired-end"){
      system ("${bwa_d}bwa mem -t $threads $human_genome ${input1}_1${suffix1} ${input1}_2${suffix1} >${input1}_h.sam");
    } elsif ($sequencing_type eq "single-end") {
      system ("${bwa_d}bwa mem -t $threads $human_genome ${input1}${suffix1} >${input1}_h.sam");
    }
    system ("${samtools_d}samtools view -bS -@ $thread_1 ${input1}_h.sam >${input1}_h.bam");
    system ("rm ${input1}_h.sam");
  }
}

##### Convert bam file to fastq file
sub convert_bamtofastq {
  my ($input1)=@_;
  my $suffix1="";
  if(($file_suffix =~ "bam" || $file_suffix =~ "sam" ) && $order ne 2){
    $suffix1=$file_suffix;
  } else {
    $suffix1="_h.bam";
  }
  if($sequencing_type eq "paired-end"){
    if(-e ${input1}.$suffix1) {
      system ("${samtools_d}samtools view -b -f 4 -F 264 -@ $thread_1 ${input1}$suffix1 > ${input1}_su.bam");
      system ("${samtools_d}samtools view -b -h -f 8 -F 260 -@ $thread_1 ${input1}$suffix1 > ${input1}_sm.bam");
    }
    if ($mode eq "standard"){
      if (-e ${input1}.$suffix1){
        system ("${samtools_d}samtools view -b -f 12 -F 256 -@ $thread_1 ${input1}$suffix1 > ${input1}_pe.bam");
 #       system ("${samtools_d}samtools view -b -f 2 -@ $thread_1 ${input1}$suffix1 > ${input1}_mpe.bam");
      }
      system("${SE_MEI_d}extractSoftclipped -l 20 ${input1}$suffix1 >${input1}_soft.fastq.gz");            ##### minumum length of split reads is set as 20 bp;
#      system("rm ${input1}_mpe.bam");
      system("${samtools_d}samtools merge -f ${input1}_h1.bam ${input1}_sm.bam ${input1}_su.bam ${input1}_pe.bam");
    } elsif($mode eq "fast"){
      system("${samtools_d}samtools merge -f ${input1}_h1.bam ${input1}_sm.bam ${input1}_su.bam");
    }
    system ("${samtools_d}samtools sort -n -@ $thread_1 ${input1}_h1.bam -o ${input1}_h1.sort.bam");
    system ("${samtools_d}samtools fastq -@ $thread_1 -N ${input1}_h1.sort.bam -1 ${input1}_h1_1.1fq -2 ${input1}_h1_2.1fq");
    system ("perl ${directory}Scripts/Check_paired_end.pl -s ${input1}_h1 -f .1fq");
    system ("mv ${input1}_h1_1.1fq2 ${input1}_h1_1.1fq");
    system ("mv ${input1}_h1_2.1fq2 ${input1}_h1_2.1fq");
    system("rm ${input1}_h1.bam");
  } elsif ($sequencing_type eq "single-end"){
    if(-e ${input1}.$suffix1){
      if ($suffix1 =~ "bam") {
        system("${samtools_d}samtools view -b -f 4 -@ $thread_1 ${input1}$suffix1 >${input1}_u.bam");
      } elsif ($suffix1 =~ "sam") {
        system("${samtools_d}samtools view -b -@ $thread_1 ${input1}_h.sam >${input1}_h.bam");
      }
    }
    system("${SE_MEI_d}extractSoftclipped -l 20 ${input1}$suffix1 >${input1}_soft.fastq.gz");              
    system("${hydra_d}bin/bamToFastq -bam ${input1}_u.bam -fq1 ${input1}_h2.1fq -fq2 ${input1}_h1.1fq");
    system ("rm ${input1}_h2.1fq");
  }
}

##### validate function
sub v_obtain_seq {
  open VIRUS_f2, "${input_sampleID}.virus_f2";
  open VIRUS_f3, ">${input_sampleID2}.virus_f2";
  while (my $tmp_v2=<VIRUS_f2>) {
    my @tmp_v2 =split /\s+/, $tmp_v2;
    my $tmp_v3 = $tmp_v2[15]."_".$tmp_v2[31]."_".$tmp_v2[32]."_".$tmp_v2[33]."_".$tmp_v2[34]."_".$tmp_v2[35];
    if ($input_sampleID2 eq $tmp_v3) {
      print VIRUS_f3 "$tmp_v2";
      last;
    }
  }
  close VIRUS_f2;
  close VIRUS_f3;
  system ("perl ${directory}Scripts/Extract_specific_loci_final_reads.pl ${input_sampleID}_f2 ${input_sampleID2}.virus_f2 >${input_sampleID2}.information");
  my $cmd=q(awk '{print$8}');
  system ("$cmd ${input_sampleID2}.information |sort |uniq >${input_sampleID2}.id");
  system ("perl ${directory}Scripts/Extract_fastq.pl -f ${input_sampleID}_1.1fuq -b ${input_sampleID2}.id -o ${input_sampleID2}_aligned_1.fuq");
  system ("perl ${directory}Scripts/Extract_fastq.pl -f ${input_sampleID}_2.1fuq -b ${input_sampleID2}.id -o ${input_sampleID2}_aligned_2.fuq");
  system ("perl ${directory}Scripts/Extract_fuq_split.pl -f ${input_sampleID}_1sf.fuq -b ${input_sampleID2}.id -o ${input_sampleID2}_aligned_sf.fuq");
  system ("perl ${directory}Scripts/Convert_fastq_to_fasta.pl ${input_sampleID2}");
  system ("perl ${directory}Scripts/Convert_fastq_to_fasta_for_split_read.pl ${input_sampleID2}");
  system ("perl ${directory}Scripts/Convert_fastq_to_fasta.pl ${input_sampleID2}");
  system ("perl ${directory}Scripts/Convert_fastq_to_fasta_for_split_read.pl $input_sampleID2");
  system ("cat ${input_sampleID2}_aligned.fas ${input_sampleID2}_aligned_sf.fas >${input_sampleID2}_aligned_both.fas");
}

sub v_bwa_validate {
  for (my $i=0;$i<@groups;$i++){
    if ($groups[$i] eq "target") {
      $ref=${directory}."Database/GI/".${GI}.".fa";
    } elsif ($groups[$i] eq "human") {
      $ref=$human_genome;
    } elsif ($groups[$i] eq "vicaller") {
      $ref=$virus_genome;
    }
    system ("${bwa_d}bwa mem -t $threads -k 19 -r 1.5 -c 100000 -m 50 -a -T 20 -h 10000 -Y -M $ref ${input_sampleID2}_aligned_1.fuq ${input_sampleID2}_aligned_2.fuq >${input_sampleID2}_aligned_c_$groups[$i]_notdefault.sam");
    system ("${samtools_d}samtools view -h -F 4 ${input_sampleID2}_aligned_c_$groups[$i]_notdefault.sam >${input_sampleID2}_aligned_c_$groups[$i]_notdefault2.sam");
    system ("rm ${input_sampleID2}_aligned_c_$groups[$i]_notdefault.sam");
    system ("${bwa_d}bwa mem -T 20 -a -t $threads $ref ${input_sampleID2}_aligned_1.fuq ${input_sampleID2}_aligned_2.fuq >${input_sampleID2}_aligned_c_$groups[$i]_default.sam");
    system ("${samtools_d}samtools view -h -F 4 ${input_sampleID2}_aligned_c_$groups[$i]_default.sam >${input_sampleID2}_aligned_c_$groups[$i]_default2.sam");
    system ("rm ${input_sampleID2}_aligned_c_$groups[$i]_default.sam");
    system ("${bwa_d}bwa mem -t $threads -k 19 -r 1.5 -c 100000 -m 50 -a -T 20 -h 10000 -Y -M $ref ${input_sampleID2}_aligned_sf.fuq >${input_sampleID2}_aligned_sf_$groups[$i]_notdefault.sam");
    system ("${samtools_d}samtools view -h -F 4 ${input_sampleID2}_aligned_sf_$groups[$i]_notdefault.sam > ${input_sampleID2}_aligned_sf_$groups[$i]_notdefault2.sam");
    system ("rm ${input_sampleID2}_aligned_sf_$groups[$i]_notdefault.sam");
    system ("${bwa_d}bwa mem -T 20 -a -t $threads $ref ${input_sampleID2}_aligned_sf.fuq >${input_sampleID2}_aligned_sf_$groups[$i]_default.sam");
    system ("${samtools_d}samtools view -h -F 4 ${input_sampleID2}_aligned_sf_$groups[$i]_default.sam >${input_sampleID2}_aligned_sf_$groups[$i]_default2.sam");
    system ("rm ${input_sampleID2}_aligned_sf_$groups[$i]_default.sam");
  }
}

sub v_combine {
  system ("cat ${input_sampleID2}_aligned_c_target_notdefault2.sam ${input_sampleID2}_aligned_sf_target_notdefault2.sam |grep -v ^@ >${input_sampleID2}_aligned_target_notdefault2.sam");
  system ("cat ${input_sampleID2}_aligned_c_target_default2.sam ${input_sampleID2}_aligned_sf_target_default2.sam |grep -v ^@ >${input_sampleID2}_aligned_target_default2.sam");
  system ("cat ${input_sampleID2}_aligned_c_vicaller_notdefault2.sam ${input_sampleID2}_aligned_sf_vicaller_notdefault2.sam |grep -v ^@ >${input_sampleID2}_aligned_vicaller_notdefault2.sam");
  system ("cat ${input_sampleID2}_aligned_c_vicaller_default2.sam ${input_sampleID2}_aligned_sf_vicaller_default2.sam |grep -v ^@ >${input_sampleID2}_aligned_vicaller_default2.sam");
  system ("cat ${input_sampleID2}_aligned_c_human_notdefault2.sam ${input_sampleID2}_aligned_sf_human_notdefault2.sam |grep -v ^@ >${input_sampleID2}_aligned_human_notdefault2.sam");
  system ("cat ${input_sampleID2}_aligned_c_human_default2.sam ${input_sampleID2}_aligned_sf_human_default2.sam |grep -v ^@ >${input_sampleID2}_aligned_human_default2.sam");
}

sub v_blat_validate {
  system ("${blat_d}blat -minIdentity=0 -minScore=0 -stepSize=4 ${directory}Database/GI/${GI}.fa ${input_sampleID2}_aligned_both.fas ${input_sampleID2}_aligned_target.psl");
  system ("${blat_d}blat -minIdentity=0 -minScore=0 -stepSize=4 $virus_genome ${input_sampleID2}_aligned_both.fas ${input_sampleID2}_aligned_vicaller.psl");
  system ("${blat_d}blat -minIdentity=0 -minScore=0 -stepSize=4 $human_genome ${input_sampleID2}_aligned_both.fas ${input_sampleID2}_aligned_human.psl");
}

sub v_blastn_validate {
  system ("${blastn_d}blastn -query ${input_sampleID2}_aligned_both.fas -word_size 10 -evalue 1e-5 -outfmt 6 -db ${directory}Database/GI/${GI}.fa >${input_sampleID2}_aligned_target.out");
  system ("${blastn_d}blastn -query ${input_sampleID2}_aligned_both.fas -word_size 10 -evalue 1e-5 -outfmt 6 -db $virus_genome >${input_sampleID2}_aligned_vicaller.out");
  system ("${blastn_d}blastn -query ${input_sampleID2}_aligned_both.fas -word_size 10 -evalue 1e-5 -outfmt 6 -db $human_genome >${input_sampleID2}_aligned_human.out");
}

sub v_convert{
  for (my $i=0;$i<@groups;$i++){
    system ("perl ${directory}Scripts/Filtered_temp.pl ${input_sampleID2}_aligned_$groups[$i]_default2.sam aaa >${input_sampleID2}_aligned_$groups[$i]_default2_infection2");
    system ("perl ${directory}Scripts/Filtered_temp2.pl ${input_sampleID2}_aligned_$groups[$i]_default2_infection2 >${input_sampleID2}_aligned_$groups[$i]_default2_infection21");
    system ("perl ${directory}Scripts/Filtered_temp.pl ${input_sampleID2}_aligned_$groups[$i]_notdefault2.sam aaa >${input_sampleID2}_aligned_$groups[$i]_notdefault2_infection2");
    system ("perl ${directory}Scripts/Filtered_temp2.pl ${input_sampleID2}_aligned_$groups[$i]_notdefault2_infection2 >${input_sampleID2}_aligned_$groups[$i]_notdefault2_infection21");
    system ("perl ${directory}Scripts/Organize_blast_all_virus.pl ${input_sampleID2}_aligned_$groups[$i] >${input_sampleID2}_aligned_$groups[$i].out2");
    system ("perl ${directory}Scripts/Organize_blat_all_virus.pl ${input_sampleID2}_aligned_$groups[$i] >${input_sampleID2}_aligned_$groups[$i].psl2");
    my $cmd=q(awk '{$4="";print$2" "$0}');
    system ("$cmd ${input_sampleID2}_aligned_$groups[$i]_default2_infection21 >${input_sampleID2}_aligned_$groups[$i]_default2_infection22");
    system ("perl ${directory}Scripts/Replace_virus_name.pl ${input_sampleID2}_aligned_$groups[$i]_default2_infection22 $directory >${input_sampleID2}_aligned_$groups[$i]_default2_infection222");
    system ("perl ${directory}Scripts/Summary_chimeric_split_reads_final031417.pl ${input_sampleID2} ${input_sampleID2}_aligned_$groups[$i]_default2_infection222 >${input_sampleID2}_aligned_$groups[$i]_default2_infection_CS");
    system ("$cmd ${input_sampleID2}_aligned_$groups[$i]_notdefault2_infection21 >${input_sampleID2}_aligned_$groups[$i]_notdefault2_infection22");
    system ("perl ${directory}Scripts/Replace_virus_name.pl ${input_sampleID2}_aligned_$groups[$i]_notdefault2_infection22 $directory >${input_sampleID2}_aligned_$groups[$i]_notdefault2_infection222");
    system ("perl ${directory}Scripts/Summary_chimeric_split_reads_final031417.pl ${input_sampleID2} ${input_sampleID2}_aligned_$groups[$i]_notdefault2_infection222 >${input_sampleID2}_aligned_$groups[$i]_notdefault2_infection_CS");
    $cmd=q(awk '{print$2" "$0}');
    system ("$cmd ${input_sampleID2}_aligned_$groups[$i].psl2 >${input_sampleID2}_aligned_$groups[$i].psl3");
    system ("perl ${directory}Scripts/Replace_virus_name.pl ${input_sampleID2}_aligned_$groups[$i].psl3 $directory >${input_sampleID2}_aligned_$groups[$i].psl222");
    system ("perl ${directory}Scripts/Summary_chimeric_split_reads_final031417.pl ${input_sampleID2} ${input_sampleID2}_aligned_$groups[$i].psl222 >${input_sampleID2}_aligned_$groups[$i].psl_CS");
    system ("$cmd ${input_sampleID2}_aligned_$groups[$i].out2 >${input_sampleID2}_aligned_$groups[$i].out3");
    system ("perl ${directory}Scripts/Replace_virus_name.pl ${input_sampleID2}_aligned_$groups[$i].out3 $directory >${input_sampleID2}_aligned_$groups[$i].out222");
    system ("perl ${directory}Scripts/Summary_chimeric_split_reads_final_blast031417.pl ${input_sampleID2} ${input_sampleID2}_aligned_$groups[$i].out222 >${input_sampleID2}_aligned_$groups[$i].out_CS");
  }
}

sub v_summary {
  system ("perl ${directory}Scripts/Summary_all_results_for_each_read_CS.pl ${input_sampleID2} _default2_infection_CS >${input_sampleID2}.default_CS");
  system ("perl ${directory}Scripts/Summary_all_results_for_each_read_CS.pl ${input_sampleID2} _notdefault2_infection_CS >${input_sampleID2}.notdefault_CS");
  system ("perl ${directory}Scripts/Summary_all_results_for_each_read_CS.pl ${input_sampleID2} .psl_CS >${input_sampleID2}.psl_CS");
  system ("perl ${directory}Scripts/Summary_all_results_for_each_read_CS.pl ${input_sampleID2} .out_CS >${input_sampleID2}.out_CS");
  system ("perl ${directory}Scripts/Summary_by_read.pl ${input_sampleID2} >${input_sampleID2}_CS2");
  system ("perl ${directory}Scripts/Extract_specific_loci_final_visualization.pl ${input_sampleID}_f2 ${input_sampleID2}.virus_f2 >${input_sampleID2}.visualization");
  system ("paste ${input_sampleID2}.default_CS ${input_sampleID2}.notdefault_CS ${input_sampleID2}.psl_CS ${input_sampleID2}.out_CS >${input_sampleID2}.CS2");
  system ("perl ${directory}Scripts/Summary_by_read_final.pl ${input_sampleID2} >${input_sampleID2}.CS3");
  system ("rm ${input_sampleID2}_aligned*");
  system ("rm ${input_sampleID2}.CS2");
  system ("rm ${input_sampleID2}_CS2");
  system ("rm ${input_sampleID2}_error");
  system ("rm ${input_sampleID2}.id");
  system ("rm ${input_sampleID2}.information");
  system ("rm ${input_sampleID2}.*_CS"); 
}

sub v_index_GI {
  if (!(-d ${directory}."Database/GI/")){
    mkdir ${directory}."Database/GI/";
  }
  if (!(-e ${directory}."Database/GI/".${GI}.".fa")){
    system ("perl ${directory}Scripts/Extract_fasta.pl $GI $virus_genome >${directory}Database/GI/${GI}.fa");
    system ("${bwa_d}bwa index -a bwtsw ${directory}Database/GI/${GI}.fa");
    system ("${blastn_d}/bin/makeblastdb -in ${directory}Database/GI/${GI}.fa -dbtype nucl");
  }
}

##### check input data
sub check_input {
  if(defined($input_sampleID) && -e ${input_sampleID}."_1".${file_suffix} && -e ${input_sampleID}."_2".${file_suffix} && $sequencing_type eq "paired-end" && ($file_suffix =~ "fq" || $file_suffix =~ "fastq")){
    print "~~~~~ paired-end reads in fastq format were loaded\n";
  } elsif (defined($input_sampleID) && -e ${input_sampleID}.${file_suffix} && $sequencing_type eq "single-end" && ($file_suffix =~ "fq" || $file_suffix =~ "fastq")){
    print "~~~~~ single-end read in fastq format was loaded\n";
  } elsif(defined($input_sampleID) && -e ${input_sampleID}.${file_suffix} && $sequencing_type eq "paired-end" && $file_suffix =~ "bam"){
    print "~~~~~ paired-end reads in bam format were loaded\n";
  } elsif (defined($input_sampleID) && -e ${input_sampleID}.${file_suffix} && $sequencing_type eq "single-end" && $file_suffix =~ "bam"){
    print "~~~~~ single-end reads in bam format were loaded\n";
  } else {
    prtErr("# Error: cound not find the input data under the provided sampleID\n");
    prtUsa();
    exit;
  }
}

##### print help page;
sub prtUsa {
  print "\n$0 <functions> [arguments]\n\n";
  print "\n";
  print "Functions:\n\tconfig|detect|validate|calculate\n\n";

  print "# Arguments for detect function:\n";
  print "	-i|input_sampleID       sample ID (required)\n";
  print "	-f|file_suffix          the suffix of the input data, including: .fq.gz|fastq.gz,.fq|fastq and .bam, indicate fastq and bam formate searately (default: .fq.gz) (required)\n";
  print "	-m|mode                 running mode, including: standard, fast (default: standard)\n";
  print "	-d|data_type		data type, including: WGS, RNA-seq (default: WGS)\n";
  print "	-s|sequencing_type	type of sequencing data, including: paired-end, single-end (default: paired-end)\n";
  print "	-t|threads		the number of threads will be used (default: 1)\n";
  print "	-r|repeat		check repeat sequence\n";
  print "	-a|align_back_to_human	reciprocal align back to the human reference genome\n";
  print "	-q|QS_cutoff		quality score cutoff for each nucleotide\n";
  print "	-c|config		user defined configure file\n";
  print "	-b|build		build version, including: hg19 and hg38 (default: hg38)\n";
  print "	-h|help			print this help\n";

  print "\n\n# Arguments for validate function:\n";
  print "	-i|input_sampleID       sample ID (required)\n";
  print "	-c|config               user defined configure file\n";
  print "	-t|threads              the number of threads will be used (default: 1)\n";
  print "	-S|String		string with sample ID, integration region, candidate virus, GI (required)\n";
  print "	-G|GI			GI (required)\n";
  print "	-V|Virus		candidate virus (required)\n";
  print "	-h|help                 print this help\n";

  print "\n\n# Arguments for calculate function:\n";
  print "	-i|input_sampleID	sample ID (required)\n";
  print "	-c|config		user defined configure file\n";
  print "	-t|threads		the number of threads will be used (default: 1)\n";
  print "	-F|File_suffix_bam	the suffix of the input bam file (required)\n";
  print "	-I|Index_sort		if the input is a sorted index BAM file\n";
  print "	-C|Chr			chromosome ID (required)\n";
  print "	-P|Position		integration site (required)\n";
  print "	-B|Breakpoint		both or one of upstream and downstream breakpoints detected, including: 1, 2 (default: 2)\n";
  print "	-N|Number_reads		number of chimeric and split reads\n";
  print "	-h|help			print the help\n";
  print "\n\n\n";

  print "Example for installing VIcaller:\n";
  print "	perl VIcaller.pl config\n\n";
  print "Examples for detecting viral integrations:\n";
  print " # WGS data in single-end fastq format:\n";
  print "	perl VIcaller.pl detect -d WGS -i seq -f .fastq.gz -m standard -s single-end -t 12 -Q 20 -a -r\n";
  print " # RNA data in paired-end fastq format:\n";
  print "	perl VIcaller.pl detect -d RNA-seq -i seq -f .fastq.gz -s paired-end -t 12\n";
  print " # RNA alignment data in bam format\n";
  print "	perl VIcaller.pl detect -d RNA-seq -i seq -f .bam -s paired-end -t 12\n\n";
  print "Examples for validating each candidate viral integration:\n";
  print "	perl VIcaller.pl validate -i seq -S seq_1_24020575_24020787_human_papillomavirus_type_218931404 -G 218931404 -V human_papillomavirus_type\n\n";
  print "Examples for calculating integration allele fraction:\n";
  print "	perl VICaller.pl calculate -i seq -f .bam -S -C 1 -P 24020575 -B 2 -N 20\n\n";
}

##### print error message with the input data;
sub prtErr{
  print "\n\n=====================================\n@_\n=====================================\n\n";
}

##### load the created config file;
sub load_config {
  if(-e ${directory}."VIcaller.config"){
    print "~~~~~ load customized config file\n";
    open CONFIG, "${directory}VIcaller.config";
    while(<CONFIG>){system("$_");}
  } else {
    prtErr("Error: could not find the VIcaller.config file under the software directory, please run the config step\n");
    exit;
  }
}



