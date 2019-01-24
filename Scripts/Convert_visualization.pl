#!usr/bin/perl
use strict;

my $line="";
my @line=();
my $read_number=0;
my %reads=();
my $read_order=0;

open F1,"$ARGV[0]";
while(<F1>){
 @line=split;
 if($line[0] eq "O1"){print"$_"; $read_number=0;}
 elsif($line[0] eq "O2"){print "$_";$read_number++;}
 elsif($line[0] eq "O3" && $line[1] =~"Start"){
  print"$_"; 
  while(<F1>){
   my @line2=split;
   if($_ =~ "human part"){print "$_";$read_order=0;%reads=();}
   elsif($_ =~ "virus part"){
    my @read_names=sort{$reads{$a}->[3]<=>$reads{$b}->[3]} keys %reads;
    for(my $i=0;$i<@read_names;$i++){print join("\t",@{$reads{$read_names[$i]}}[0..2]);print "\n";}
    $read_order=0;%reads=();
    print "$_";
                           }
   elsif($line2[0] =~ "Seq"){
    if(exists($reads{$line2[0]})){${$reads{$line2[0]}}[1]=${$reads{$line2[0]}}[1].$line2[2];${$reads{$line2[0]}}[2]=$line2[3];}
    unless(exists($reads{$line2[0]})){@{$reads{$line2[0]}}=@line2[1..3];push(@{$reads{$line2[0]}},$read_order);$read_order++;}
                            }
   elsif($line2[0] eq "O3" && $line2[1] =~ "End"){
    my @read_names=sort{$reads{$a}->[3]<=>$reads{$b}->[3]} keys %reads;
    for(my $i=0;$i<@read_names;$i++){print join("\t",@{$reads{$read_names[$i]}}[0..2]);print "\n";}
    $read_order=0;%reads=();
    print "$_";
    last;
                                                 }
             }
                                              }
           }
