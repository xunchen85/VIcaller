#!usr/bin/perl
use strict;
my $line="";
my @line=();
my %mapped=();
my $distance=0;
my $direction="";
my %single_mapped=();
my $name="";
my $name2="";
my $start1=0;
my $start2=0;
my $length=0;
my $end=0;
my @line1=();
my @line2=();

open FILE,"$ARGV[0]_$ARGV[1].sam";
open MAPPED,">$ARGV[0]_$ARGV[1].hmap";
open UNMAPPED,">$ARGV[0]_$ARGV[1].hunmapped";
while(<FILE>){
 @line=split;
 $line[13]=~s/AS:i://;
 $line[12]=~s/MD:Z://;

 $start1=0;$end=0;$length=0;
 if($line[0]=~"@"){next;}
 my @temp1=($line[5]=~/(\d+)/g);
 my @temp2=($line[5]=~/([A-Z])/g);
 for(my $i=0;$i<@temp1;$i++){
  if($temp2[$i] eq "S"){
   if($i==0){$start1=$temp1[$i]+1;$end=$temp1[$i];$length=$temp1[$i];}
   if($i==@temp1-1){$length+=$temp1[$i];}
                       }
  elsif($temp2[$i] eq "D"){next;}
  elsif($temp2[$i] eq "M" || $temp2[$i] eq "I"){
   if($i==0){$start1=1;$end=$temp1[$i];$length=$temp1[$i];}
   elsif($i==@temp1-1){$length+=$temp1[$i];$end+=$temp1[$i];}
   else{$end+=$temp1[$i];$length+=$temp1[$i]}  }
                            }
#####
 if($line[1]%256>=128){$direction="R";}
 elsif($line[1]%256<128 && $line[1]%256>=64){$direction="L";}

#####
if($line[1]%8>=4){
  print UNMAPPED "$direction $line[0] 0 $start1 $end $length 0 0 0 $direction @line[1..3] @line[5..8] $line[13]\n";
                 }
elsif($line[1]%16>=8){
  if($direction eq "L"){
    print MAPPED "$direction $line[0] 0 $start1 $end $length 0 0 0 $direction @line[1..3] @line[5..8] $line[13]\n";
                       }
  else{
    print MAPPED "$direction $line[0] 0 0 0 0 $start1 $end $length $direction @line[1..3] @line[5..8] $line[13]\n";
      }
                     }
#####
else{
  $name=$line[0]."_".$direction."_".$line[2]."_".$line[3];
  if($line[6] eq "="){$line[6]=$line[2];}
  $mapped{$name}=join(" ",($line[0],$start1,$end,$length,$direction,@line[1..3],@line[5..8],$line[13]));
    }
              }
#########
my @name_mapped=keys %mapped; my $note="";
  @line2=(); @line1=();
for(my $i=0;$i<@name_mapped;$i++){
 @line1=split /\s+/, $mapped{$name_mapped[$i]};
 if($line1[4] eq "L"){
  $name2=$line1[0]."_R_".$line1[9]."_".$line1[10];
  if(exists($mapped{$name2})){                           ##### calculate the distance and change note; colume 2 and 3;
    @line2=split /\s+/, $mapped{$name2};
    if($line1[6] eq $line2[6]){
     if($line1[10] eq $line2[7] && $line1[7] eq $line2[10] && abs($line1[11]) eq abs($line2[11])){
       $note="P1";$distance=abs($line1[11]);}
     else{
      $note="P2";
      my $te1=0;
      if($line1[7]<$line2[7]){$te1=$line1[3];}
      else{$te1=$line2[3];}
      
      $distance=abs($line2[7]-$line1[7])+$te1;
         }
                              }
    else{$distance=0; $note="P3";}

#############
    print MAPPED "$note $line1[0] $distance @line1[1..3] @line2[1..3] @line1[4..12] @line2[4..12]\n";
                             }
  else{
   print MAPPED "E $line1[0] 0 0 0 0 0 0 0 @line1[4..12]\n";}
                    }
 else{
  $name2=$line1[0]."_L_".$line1[9]."_".$line1[10];
  if(exists($mapped{$name2})){
    @line2=split /\s+/, $mapped{$name2};
    if($line1[6] eq $line2[6]){
     if($line1[10] eq $line2[7] && $line1[7] eq $line2[10] && abs($line1[11]) eq abs($line2[11])){$note="P1";$distance=abs($line1[11]);}
     else{
      $note="P2";
      my $te1=0;
      if($line1[7]<$line2[7]){$te1=$line1[3];}
      else{$te1=$line2[3];}
      $distance=abs($line2[7]-$line1[7])+$te1;
         }
                                         }
    else{
     $distance=0; $note="P3";
        }

#############
   print MAPPED "$note $line1[0] $distance @line2[1..3] @line1[1..3] @line2[4..12] @line1[4..12]\n";
                             }
  else{
   print MAPPED "E $line1[0] 0 0 0 0 0 0 0 @line1[4..12]\n";}
      }
                                 }
close MAPPED;
close UNMAPPED;
