#!usr/bin/perl
use strict;

my %list=();
my @list=();
my @line=();
my $line="";
my $title="";
my $title2="";
my $sample=$ARGV[0];
my $ok=0;

##########################
open F1,"$ARGV[1]";
while(<F1>){
 @line=split;
 $line[0]=~s/ad://;
 $title=$line[0]."_".$line[15]."_".$line[35]."_".$line[30]."_".$line[31]."_".$line[32]."_".$line[33];
 $list{$title}=$line[34];
}

##########################
open F2,"$ARGV[0]";
my @title_O1=();
my $title3="";
while(<F2>){
 @line=split;
  if($line[0] eq "O1"){
   $ok=0;
   $title2=$line[2]."_".$line[1]."_".$line[22]."_".$line[17]."_".$line[18]."_".$line[19]."_".$line[20];
   $title3="";
   $line[21]= lc $line[21];
   if($line[21] =~ "unknown"){$line[21]=$line[22];}
   if(exists($list{$title2})){
    $ok=1;
    @title_O1=();
    @title_O1=@line;
    $title3=$line[1]." ".$line[18]." ".$line[19]." ".$line[20]." ".$list{$title2}." ".$line[22];
   }
                      }
   elsif($line[0] eq "O2" && $ok==1){
   my $title4=$title3." ".$line[0]." ".$line[1]." ".$line[2]." ".$line[3]." ".$line[4]." ".$line[5]." ".$line[7]." ".$line[9];
   print "$title4\n";
  }
           }
