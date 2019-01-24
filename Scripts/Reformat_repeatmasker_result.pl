#!usr/bin/perl
use strict;
my $line="";
my @line=();
my %read=();

open REPEAT,"$ARGV[0].fasta3.out";
while(<REPEAT>){
 @line=split;
 unless($line[0]){next;}
 if($line[0]=~"[a-zA-Z]"){next;}
 if(exists($read{$line[4]})){
  push(@{$read{$line[4]}},@line[5..6]);push(@{$read{$line[4]}},@line[9..10]);
                            }
 else{
  push(@{$read{$line[4]}},@line[5..6]);push(@{$read{$line[4]}},@line[9..10]);
     }
               }
@line=keys %read;
print "Read_name Start End Matching_repeat Repeat_class_family\n";
for(my $i=0;$i<@line;$i++){
 print"$line[$i] @{$read{$line[$i]}}\n";
}
