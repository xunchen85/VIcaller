#!usr/bin/perl
use strict;

my $line="";
my @line=();
my $type1="";
my $k=0;
my @len=();
my $as="";

###############
open F1,"$ARGV[0]";

while(<F1>){
 $as="";
 my @line=split;
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

###############
 for(my $i=11;$i<@line;$i++){
  if($line[$i]=~"AS"){
   $as=$line[$i];
   $as=~s/AS\:i\://;
                     }
                            }
###############
                 }
 else{$type1="U";next;}
 my $virus_name="";
  print "$type1 @line[0..8] $as $as\n";
}
