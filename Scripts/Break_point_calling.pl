#!usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;

################
my $position="test_sm.sam";
my $virus="test_vsu.sam";
my $type="test.type";
my $output="test";
my $taxonomy="virus_db.taxonomy";

GetOptions('position=s'=>\$position,
           'virus=s'=>\$virus,
           'type=s'=>\$type,
           'taxonomy=s'=>\$taxonomy,
           'output=s'=>\$output);

################
my %position=();
my @position=();
my @virus=();
my @type=();
my %type=();
my $md="";
my $as="";

#################
open TYPE, "$type";
while(<TYPE>){
 @type=split;
 if($type[1] eq "L" || $type[1] eq "R"){$type{$type[0]}=$type[1];}
             }

#################
open TAXONOMY, "$taxonomy";

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

#################    # Read position information;
open POSITION, "$position";
while(<POSITION>){
 @position=split;
 $md="";$as="";
 for(my $i=11;$i<@position;$i++){
  if($position[$i]=~s/MD:Z://){$md=$position[$i];}
  if($position[$i]=~s/AS:i://){$as=$position[$i];}
                                }
 $position{$position[0]}=$position[0]." ".$position[1]." ".$position[2]." ".$position[3]." ".$position[4]." ".$position[5]." ".$position[6]." ".$position[7]." ".$position[8]." ".$md." ".$as;
                 }
##################
open VIRUS, "$virus";
open OUT1, ">${output}_breakpoint";
while(<VIRUS>){
 @virus=split;
 $md="NA";$as="NA";
 if($virus[2] eq "*"){next;}
 for(my $i=11;$i<@virus;$i++){
  if($virus[$i]=~s/MD:Z://){$md=$virus[$i];}
  if($virus[$i]=~s/AS:i://){$as=$virus[$i];}
                             }
 if($as<22){next;}
 if(exists($type{$virus[0]})){
  if(exists($taxonomy{$virus[2]})){print OUT1 "$type{$virus[0]} $taxonomy{$virus[2]}[10] $position{$virus[0]} @virus[1..8] $md $as\n";}
  else{print OUT1 "$type{$virus[0]} unknown $position{$virus[0]} @virus[1..8] $md $as\n";}
                             }
              }
