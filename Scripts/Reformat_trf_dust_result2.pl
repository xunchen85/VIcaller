#!/usr/bin/perl
use strict;
use warnings;

####################### read repeatmasker repeat file
open REPEAT, "$ARGV[0].repeat";
my %repeats=();
my @repeat_masker=();
my $line2=<REPEAT>;

print "$line2";                            # title line
while(<REPEAT>){
 my @line=split;
 chomp @line;
 @{$repeats{$line[0]}}=@line;
               }
###################### read file of trf/dust results, make it hash;
open FASTA, "$ARGV[0].fasta2";
my %fasta=();
my $name="";
while(<FASTA>){
 if($_=~"^>"){
  $name=$_;
  $name=~s/^>//;
  chomp($name);
  $fasta{$name}="";
             }
  else{
   my $seq=$_;
   chomp($seq);
   $fasta{$name}=$fasta{$name}.$seq;
      }
               }

###################### read fasta2 hash
my @seq=();
my @fasta=();
my @position=();
my $region=0;
my @list=keys %fasta;

###########
for(my $d=0;$d<@list;$d++){      # loop to read each read of trf/dust results start
  @seq=();
  @fasta=();
  @position=();
  $region=0;
  my $line=$list[$d];
  @seq=split //, $fasta{$line};
  unless($fasta{$line}=~/N/i){next;}
  @repeat_masker=();
  if(exists($repeats{$line})){                      # if repeatmasker already detect this seq contained repeats;
   @repeat_masker=@{$repeats{$line}};               # repeat masker results;

###########
   for(my $i=0;$i<@seq;$i++){                       # read each nr of this read;
    if($seq[$i] eq "N" || $seq[$i] eq "n") {
     my $not_in=0;
     for(my $j=1;$j<@repeat_masker;$j+=4){
      if($i+1<=$repeat_masker[$j+1] && $i+1>=$repeat_masker[$j]){$j=$#repeat_masker+2;$not_in=-1;}
                                         }
     if($not_in==-1){next;}
     else{
      $region=0;
      if(@position){
        if($i+1<$position[$#position-2]+1){next;}
        elsif($i+1==$position[$#position-2]+1){$position[$#position-2]=$i+1;}
        elsif($i+1>$position[$#position-2]+1){push(@position,($i+1,$i+1,"trf","dust"));}
                   }
      else{push(@position,($i+1,$i+1,"trf","dust"));}
         }
                                           }
                            }
########### merge two repeatmasker and trf/dust results
 for(my $i=0;$i<@position;$i+=4){
  if($position[$i+1]<$repeat_masker[1]){
    splice @repeat_masker,1,0,@position[$i..$i+3];
                                       }
  elsif($position[$i]>$repeat_masker[$#repeat_masker-2]){push(@repeat_masker,@position[$i..$i+3]);}
  else{
   for(my $j=1;$j<@repeat_masker;$j+=4){
                                       }
      }
                                 }
   @{$repeats{$line}}=@repeat_masker;               # combined all repeats;
                             }                      # end of repeatmasker results
  else{                                             # if repeatmasker did not detect this seq contained repeats;
   for(my $i=0;$i<@seq;$i++){                      # read each nr of this read;
    if($seq[$i] eq "N" || $seq[$i] eq "n") {
      if(@position){
        if($i+1<$position[$#position-2]+1){next;}
        elsif($i+1==$position[$#position-2]+1){$position[$#position-2]=$i+1;}
        elsif($i+1>$position[$#position-2]+1){push(@position,($i+1,$i+1,"trf","dust"));}
                   }
      else{push(@position,$i+1,$i+1,"trf","dust");}

                                           }
                            }
   @{$repeats{$line}}=($line,@position);            # only trf/dust repeats;
      }                                             # end of not repeatmasker
}

#################### output repeat2 file
my @list2=keys %repeats;
for(my $i=0;$i<@list2;$i++){
 chomp @{$repeats{$list2[$i]}};
 print "@{$repeats{$list2[$i]}}\n";
}


