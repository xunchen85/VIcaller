########################## extract clean paired-end unmapped reads.
#!usr/bin/perl
use strict;
my $directory=$ARGV[1];

open S3,">$ARGV[0]_aligned_sf.fas";
open F1,"$ARGV[0]_aligned_sf.fuq"||di("di too");
my @line=split;
my $read_name="";

while($read_name=<F1>){
   my $read_1=$read_name;
   $read_1=~s/^@/>/g;
   my @temp11=split /\s+/, $read_1;
   $read_name=$temp11[0];
   $read_name=~s/\/1//;
   $read_name=~s/^@//;
   my $temp01_1=<F1>;
   my $temp02_1=<F1>;
   my $temp03_1=<F1>;
   my $temp04_1=$read_1.$temp01_1;
   print S3 "$temp04_1";
                      }
