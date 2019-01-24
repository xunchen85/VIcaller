#!usr/bin/perl
use strict;
my $line="";
my @line=();
my $type1="";
my $k=0;
my @len=();
my $as="";
my $as_min=$ARGV[2];

##############
open F2,"$ARGV[1]";                            ##### type file to remove chimeric and split reads;
my %type=();
while(<F2>){
 @line=split;
 if($line[1] eq "L" || $line[1] eq "R"){$type{$line[0]}=$line[1];}
}

##############
open F1,"$ARGV[0]";
while(<F1>){
 $as="";
 my @line=split;
 if(exists($type{$line[0]})){next;}
 if($line[0]=~"@"){next;}
 @len=($line[5]=~/(\d+)/g);
 my $len=0; map {$len+=$_} @len;
 $k=0;my $type1="";
 if($line[1]%8<4){ 
 if($line[1]%4>=2 && $line[1]%256>=128){$type1="PR";$k=1;}
 elsif($line[1]%4>=2 && $line[1]%256<128){$type1="PF";$k=1;}
 elsif($line[1]%256<128){$type1="F";$k=1;}
 elsif($line[1]%256>=128){$type1="R";$k=1;}
 else{$type1="N";$k=1;}

##############
 for(my $i=11;$i<@line;$i++){
  if($line[$i]=~"AS"){
   $as=$line[$i];
   $as=~s/AS\:i\://;
                     }
                            }
##############
                 }
 else{$type1="U";next;}
 if($as<$as_min){next;}                             ##### For Chimeric and split reads, we will keep those results

##############
 my $virus_name="";
  print "$type1 @line[0..8] $as $as\n";
}
