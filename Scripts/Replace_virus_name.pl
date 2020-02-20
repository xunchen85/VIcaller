#!/usr/bin/perl
use strict;
############################## This script did not exclude some unknow viral sequences.
##### modify to change all names related to HBV

my $line="";
my @line=();
my %taxo=();
my $virus_list=$ARGV[1];
open F2,"$virus_list";
while(<F2>){
 @line=split;
 $taxo{$line[0]}=$line[1];
}


open F1,"$ARGV[0]";
while(<F1>){
 @line=split;
 if(exists($taxo{$line[0]})){
  if($taxo{$line[0]} =~ "hbv" || $taxo{$line[0]} =~ "hepatitis_b"){$line[0]="hepatitis_b_virus";}
  else{$line[0]=$taxo{$line[0]};}

#  $line[0]=$taxo{$line[0]};
  print "@line\n";
                            }
 elsif(!(exists($taxo{$line[0]})) && !($ARGV[0] =~ "human")){
  print "@line\n";
 }
 elsif($ARGV[0] =~ "human"){
  if($line[0] eq "NC_007605"){next;}
  print "@line\n";
 }
}
