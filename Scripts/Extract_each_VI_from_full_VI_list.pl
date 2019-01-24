#!use/bin/perl
use strict;
my $line="";
my @line=();
my $i=0;
my $cpu=$ARGV[2];
my $Main_folder=$ARGV[1];

########################
my %taxo=();
open DB,"${Main_folder}Scripts/virus_list_for_virus_c_contamination_change_HHV4_type_2_to_HHV4";
while(<DB>){
 @line=split;
 $taxo{$line[0]}=$line[1];
}

########################
open FILE,"$ARGV[0]";
while(<FILE>){
 @line=split;
 my $file_name="";
 my $virus="N";
 if(exists($taxo{$line[35]})){
  if($taxo{$line[35]} =~ "hbv" || $taxo{$line[35]} =~ "hepatitis_b"){$virus="hepatitis_b_virus";}
  else{$virus=$taxo{$line[35]};}
  $line[34]=$virus;
  $file_name=$line[15]."_".$line[31]."_".$line[32]."_".$line[33]."_".$line[35]."_".$virus;
 }
 else{
  $file_name=$line[15]."_".$line[31]."_".$line[32]."_".$line[33]."_".$line[34]."_".$line[35];
 }
 open OUTPUT,">$file_name.virus_f2";
 print OUTPUT "@line[0..$#line]\n";
 close OUTPUT;
 print "./Checking_integration_read_level_final081817-2.sh $line[15] $file_name $line[35] $virus $cpu\n";
             }
