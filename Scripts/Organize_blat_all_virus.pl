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
open BLAT,"$ARGV[0].psl";

$line=<BLAT>;
$line=<BLAT>;$line=<BLAT>;$line=<BLAT>;$line=<BLAT>;
while(<BLAT>){                                      ### read each line
 @line=split;
 $read_title=$line[9];
 $read_title=~s/\/[12]$//;
 unless($array_title){$array_title=$read_title;}
 if($array_title && $array_title eq $read_title || $#array <0 ){    ### if the same array (the same paired-end read)
  push(@array,[@line]);
  if(eof(BLAT)){$array_title=$read_title; goto LAST;}
                                                }
 elsif(!(eof(BLAT))){                                               ### first_line of each array;
  ################### Organize each read
LAST:  
  for(my $i=0;$i<@array;$i++){
   if($array[$i][9] =~ "\/1"){push(@array_left,[@{$array[$i]}]);}
   else{
    my $paired_or_not=0;
    for(my $j=0;$j<@array_left;$j++){
     if((abs($array[$i][15]-$array_left[$j][15])<1000 && $array[$i][13] eq $array_left[$j][13]) && (($array[$i][8] eq "+" && $array_left[$j][8] eq "-") || (($array[$i][8] eq "-" && $array_left[$j][8] eq "+")))){push(@{$array_left[$j]},@{$array[$i]});$paired_or_not=1;}
                                    }
    if($paired_or_not eq 0){push(@array_right,[@{$array[$i]}]);}
       }
  }
  for(my $i=0;$i<@array_left;$i++){
   if(@{$array_left[$i]}>30){
    my $evalue=($array_left[$i][0]+$array_left[$i][21]);
    print"P $array_left[$i][13] $evalue $array_title @{$array_left[$i]}\n";}
   else{print"F $array_left[$i][13] $array_left[$i][0] $array_title @{$array_left[$i]}\n";}
                                  }
  for(my $i=0;$i<@array_right;$i++){print"R $array_right[$i][13] $array_right[$i][0] $array_title @{$array_right[$i]}\n";}
  @array_output=();@array_left=();@array_right=();
  $array_title=$read_title;
  @array=();
  push(@array,[@line]);
 }

####### eof
 elsif(eof(BLAT)){
  for(my $i=0;$i<@array;$i++){
   if($array[$i][9] =~ "\/1"){push(@array_left,[@{$array[$i]}]);}
   else{
    my $paired_or_not=0;
    for(my $j=0;$j<@array_left;$j++){
     if((abs($array[$i][15]-$array_left[$j][15])<1000 && $array[$i][13] eq $array_left[$j][13]) && (($array[$i][8] eq "+" && $array_left[$j][8] eq "-") || (($array[$i][8] eq "-" && $array_left[$j][8] eq "+")))){push(@{$array_left[$j]},@{$array[$i]});$paired_or_not=1;}
                                    }
    if($paired_or_not eq 0){push(@array_right,[@{$array[$i]}]);}
       }
  }
  for(my $i=0;$i<@array_left;$i++){
   if(@{$array_left[$i]}>30){
    my $evalue=($array_left[$i][0]+$array_left[$i][21]);
    print"P $array_left[$i][13] $evalue $array_title @{$array_left[$i]}\n";}
   else{print"F $array_left[$i][13] $array_left[$i][0] $array_title @{$array_left[$i]}\n";}
                                  }
  for(my $i=0;$i<@array_right;$i++){print"R $array_right[$i][13] $array_right[$i][0] $array_title @{$array_right[$i]}\n";}
  @array_output=();@array_left=();@array_right=();
  $array_title=$read_title;
  @array=();  
  push(@array,[@line]);
for(my $i=0;$i<@array;$i++){
   if($array[$i][9] =~ "\/1"){push(@array_left,[@{$array[$i]}]);}
   else{
    my $paired_or_not=0;
    for(my $j=0;$j<@array_left;$j++){
     if((abs($array[$i][15]-$array_left[$j][15])<1000 && $array[$i][13] eq $array_left[$j][13]) && (($array[$i][8] eq "+" && $array_left[$j][8] eq "-") || (($array[$i][8] eq "-" && $array_left[$j][8] eq "+")))){push(@{$array_left[$j]},@{$array[$i]});$paired_or_not=1;}
                                    }
    if($paired_or_not eq 0){push(@array_right,[@{$array[$i]}]);}
       }
  }
  for(my $i=0;$i<@array_left;$i++){
   if(@{$array_left[$i]}>30){
    my $evalue=($array_left[$i][0]+$array_left[$i][21]);
    print"P $array_left[$i][13] $evalue $array_title @{$array_left[$i]}\n";}
   else{print"F $array_left[$i][13] $array_left[$i][0] $array_title @{$array_left[$i]}\n";}
                                  }
  for(my $i=0;$i<@array_right;$i++){print"R $array_right[$i][13] $array_right[$i][0] $array_title @{$array_right[$i]}\n";}
}
}

