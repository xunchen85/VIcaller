#!usr/bin/perl
use strict;
my %each_virus=();
my @line=();
my $candidate_virus="";
my @top_GI=();                    ##### column 13
my @top_Unique_reads=();          ##### column 4 
my @top_Virus_length=();          ##### column 5
my $min_reads=$ARGV[1];
my $end=0;

open F1,"$ARGV[0]";
while(<F1>){
 @line=split;
 if($line[5]<$min_reads && $line[12]<$min_reads){next;}                           ##### Removed virus candidate with no more than a total number of five reads in default;
 my $temp1=$line[0]."|".$line[1];
 if($candidate_virus eq "" || $candidate_virus eq $temp1){
   $candidate_virus=$temp1;
   my $name=$line[0]."|".$line[1]."|".$line[2];
   @{$each_virus{$name}}=@line;
   if(eof(F1)){$end=1;goto LAST;}
                                                         }
 else{
LAST:

############ Top number of unique reads of the GI;
  my @name2=sort{$each_virus{$b}->[3] <=>$each_virus{$a}->[3]} keys %each_virus;
  for(my $i=0;$i<@name2;$i++){
   if(@top_Unique_reads){
    if($top_Unique_reads[3] eq ${$each_virus{$name2[$i]}}[3]){
     if($top_Unique_reads[12] < ${$each_virus{$name2[$i]}}[12]){@top_Unique_reads=@{$each_virus{$name2[$i]}};}      #### if this GI have more reads
     elsif($top_Unique_reads[12] eq ${$each_virus{$name2[$i]}}[12] && $top_Unique_reads[4] < ${$each_virus{$name2[$i]}}[4]){
      @top_Unique_reads=@{$each_virus{$name2[$i]}};}                                                                #### if this GI have the same reads, but longer genome length;  
                                                             }
                        }
   else{@top_Unique_reads=@{$each_virus{$name2[$i]}};}
                             }
############ Top number of reads of the GI;
  my @name2=sort{$each_virus{$b}->[12] <=>$each_virus{$a}->[12]} keys %each_virus;
  for(my $i=0;$i<@name2;$i++){
   if(@top_GI){
    if($top_GI[12] eq ${$each_virus{$name2[$i]}}[12]){
     if($top_GI[3] < ${$each_virus{$name2[$i]}}[3]){@top_Unique_reads=@{$each_virus{$name2[$i]}};}      #### if this GI have top unique reads
     elsif($top_GI[3] eq ${$each_virus{$name2[$i]}}[3] && $top_GI[4] < ${$each_virus{$name2[$i]}}[4]){
      @top_GI=@{$each_virus{$name2[$i]}};}                                                    #### if this GI have the same unique reads, but longer genome length;  
                                                             }
                        }
   else{@top_GI=@{$each_virus{$name2[$i]}};}
                             }
############ Top Virus length of the GI;
  my @name2=sort{$each_virus{$b}->[4] <=>$each_virus{$a}->[4]} keys %each_virus;
  for(my $i=0;$i<@name2;$i++){
   if(@top_Virus_length){
    if($top_Virus_length[4] eq ${$each_virus{$name2[$i]}}[4]){
     if($top_Virus_length[12] < ${$each_virus{$name2[$i]}}[12]){@top_Virus_length=@{$each_virus{$name2[$i]}};}      #### if this GI have top number of reads;
     elsif($top_Virus_length[12] eq ${$each_virus{$name2[$i]}}[12] && $top_Virus_length[3] < ${$each_virus{$name2[$i]}}[3]){
      @top_Virus_length=@{$each_virus{$name2[$i]}};}                                                    #### if this GI have the same number of reads, but more uniue reads;  
                                                             }
                        }
   else{@top_Virus_length=@{$each_virus{$name2[$i]}};}
                             }
  print "@top_Virus_length\n";
  print "@top_GI\n";
  print "@top_Unique_reads\n";
  @top_Virus_length=();@top_GI=();@top_Unique_reads=();

############
  %each_virus=();
  $candidate_virus=$temp1;
  my $name=$line[0]."|".$line[1]."|".$line[2];
  @{$each_virus{$name}}=@line;
  if(eof(F1) && $end eq 0){$end=1; goto LAST;}
 }
}
