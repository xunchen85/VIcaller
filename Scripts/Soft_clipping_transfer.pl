#!usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;

my $file="test_vsoft.sort.sam";
my $output="test_vsoft_breakpoint";
my $taxonomy="virus_db.taxonomy";

GetOptions('file=s'=>\$file,
           'output=s'=>\$output,
           'taxonomy=s'=>\$taxonomy
          );

open FILE, "$file";
open OUTPUT, ">$output";
open TAXONOMY, "$taxonomy";

my @line=();
my @name=();
my $md="";
my $as="";

#####################
my %taxonomy=();
sub taxon{
while(<TAXONOMY>){
 my @information=split /\$/, $_;
 my @linear=split /\;/, $information[3];
 $linear[0]=~s/ /_/g;
 $taxonomy{$information[0]}=[$information[1],$information[2],$linear[0],$linear[1],$linear[2],$linear[3],$linear[4],$linear[5],$linear[6],$information[3],$information[9]];
           }
         }
taxon();
#####################
while(<FILE>){
 @line=split;
 @name=split /\|/, $line[0];
 $as="NA";$md="NA";
 unless($line[2] eq "*"){
  for(my $i=11;$i<@line;$i++){
   if($line[$i]=~s/MD:Z://){$md=$line[$i];}
   if($line[$i]=~s/AS:i://){$as=$line[$i];}
                             }
   if($as<20){next;}
   if(exists($taxonomy{$line[2]})){
     print OUTPUT "S1 $taxonomy{$line[2]}[10] @name[1..11] @line[1..8] $md $as\n";}
   else{print OUTPUT "S1 unknown @name[1..11] @line[1..8] $md $as\n";}
                        }
              }
