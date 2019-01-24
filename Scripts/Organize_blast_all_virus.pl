#!usr/bin/perl
use strict;
my $line="";
my @line=();
my $read_title="";
my @array=();
my $array_title="";
my @array_output=();
my @array_left=();
my @array_right=();
my $top_evalue=0;    ### edit on 8/19/17
my %count_hit=0;     ### edit on 8/19/17

open BLAST,"$ARGV[0].out";
while(<BLAST>){                                      ### read each line
 @line=split;
 $read_title=$line[0];
 $read_title=~s/\/[12]$//;
 unless($array_title){$array_title=$read_title;}
 if($array_title && $array_title eq $read_title || $#array <0 ){    ### if the same array (the same paired-end read)
  if(exists($count_hit{$line[0]})){${$count_hit{$line[0]}}[0]++;}
  else{${$count_hit{$line[0]}}[0]=0;${$count_hit{$line[0]}}[1]=$line[10];}
  if($ARGV[0] =~ "_aligned_human" && ${$count_hit{$line[0]}}[0] >=1000 && ${$count_hit{$line[0]}}[1] < $line[10]){next;}
  else{push(@array,[@line]);}
  if(eof(BLAST)){$array_title=$read_title; goto LAST;}

                                                               }
 elsif(!(eof(BLAST))){                                               ### first_line of each array;

  ################### Organize each read
LAST:
  for(my $i=0;$i<@array;$i++){
   if($array[$i][0] =~ "\/1"){push(@array_left,[@{$array[$i]}]);}
   else{
    my $paired_or_not=0;
    for(my $j=0;$j<@array_left;$j++){
     if((abs($array[$i][8]-$array_left[$j][8])<1000 && $array[$i][1] eq $array_left[$j][1]) && ((($array[$i][8] - $array[$i][9] <= 0) && ($array_left[$j][8] - $array_left[$j][9] >= 0)) || (
($array[$i][8] - $array[$i][9] >= 0) && ($array_left[$j][8] - $array_left[$j][9] <= 0)))){push(@{$array_left[$j]},@{$array[$i]});$paired_or_not=1;}
                                    }
    if($paired_or_not eq 0){push(@array_right,[@{$array[$i]}]);}
       }
  }

  for(my $i=0;$i<@array_left;$i++){
   if(@{$array_left[$i]}>20){
    my $evalue=($array_left[$i][10]+$array_left[$i][22])/2;
    print"P $array_left[$i][1] $evalue $array_title @{$array_left[$i]}\n";}
   else{print"F $array_left[$i][1] $array_left[$i][10] $array_title @{$array_left[$i]}\n";}
                                  }
  for(my $i=0;$i<@array_right;$i++){print"R $array_right[$i][1] $array_right[$i][10] $array_title @{$array_right[$i]}\n";}
  @array_output=();@array_left=();@array_right=();
  $array_title=$read_title;
  @array=();
  push(@array,[@line]);
                     }
 elsif(eof(BLAST)){
  for(my $i=0;$i<@array;$i++){
   if($array[$i][0] =~ "\/1"){push(@array_left,[@{$array[$i]}]);}
   else{
    my $paired_or_not=0;
    for(my $j=0;$j<@array_left;$j++){ 
     if((abs($array[$i][8]-$array_left[$j][8])<1000 && $array[$i][1] eq $array_left[$j][1]) && ((($array[$i][8] - $array[$i][9] <= 0) && ($array_left[$j][8] - $array_left[$j][9] >= 0)) || (
($array[$i][8] - $array[$i][9] >= 0) && ($array_left[$j][8] - $array_left[$j][9] <= 0)))){push(@{$array_left[$j]},@{$array[$i]});$paired_or_not=1;}
                                    }
    if($paired_or_not eq 0){push(@array_right,[@{$array[$i]}]);}
       }
  }  
  for(my $i=0;$i<@array_left;$i++){
   if(@{$array_left[$i]}>20){
    my $evalue=($array_left[$i][10]+$array_left[$i][22])/2;
    print"P $array_left[$i][1] $evalue $array_title @{$array_left[$i]}\n";}
   else{print"F $array_left[$i][1] $array_left[$i][10] $array_title @{$array_left[$i]}\n";}
                                  }
  for(my $i=0;$i<@array_right;$i++){print"R $array_right[$i][1] $array_right[$i][10] $array_title @{$array_right[$i]}\n";}
  @array_output=();@array_left=();@array_right=();
  $array_title=$read_title;
  @array=();
  push(@array,[@line]);
  for(my $i=0;$i<@array;$i++){
   if($array[$i][0] =~ "\/1"){push(@array_left,[@{$array[$i]}]);}
   else{
    my $paired_or_not=0;
    for(my $j=0;$j<@array_left;$j++){
     if((abs($array[$i][8]-$array_left[$j][8])<1000 && $array[$i][1] eq $array_left[$j][1]) && ((($array[$i][8] - $array[$i][9] <= 0) && ($array_left[$j][8] - $array_left[$j][9] >= 0)) || (
($array[$i][8] - $array[$i][9] >= 0) && ($array_left[$j][8] - $array_left[$j][9] <= 0)))){push(@{$array_left[$j]},@{$array[$i]});$paired_or_not=1;}
                                    }
    if($paired_or_not eq 0){push(@array_right,[@{$array[$i]}]);}
       }
  }
  for(my $i=0;$i<@array_left;$i++){
   if(@{$array_left[$i]}>20){
    my $evalue=($array_left[$i][10]+$array_left[$i][22])/2;
    print"P $array_left[$i][1] $evalue $array_title @{$array_left[$i]}\n";}
   else{print"F $array_left[$i][1] $array_left[$i][10] $array_title @{$array_left[$i]}\n";}
                                  }
  for(my $i=0;$i<@array_right;$i++){print"R $array_right[$i][1] $array_right[$i][10] $array_title @{$array_right[$i]}\n";} 
              } 
}
