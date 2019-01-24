#!usr/bin/perl
use strict;

my $line="";
my @line=();
my %reads=();
my $name="";
my $name2="";
my %paired_end=();                     #### both end mapped to the same virus and in proper pair;
my %paired_end2=();                    #### both end mapped to the same virus, but not proper pair;
my %single_left=();
my %single_right=();
my %reads=();

#####################
open F1,"$ARGV[0]";
while(<F1>){           ### 1. loop start
 @line=split /\s+/, $_;
 $name=$line[4]."|".$line[8];
 unless(exists($reads{$line[4]})){$reads{$line[4]}=0;}
 if($line[0] eq "P"){        ### paired end in proper pair
  if(!(exists($paired_end{$line[4]}))){push(@{$paired_end{$line[4]}},(1,1,$line[2],"P",$line[8])); push(@{$paired_end{$name}},($line[1],$line[2],$line[7]));}
  elsif(exists($paired_end{$line[4]}) && exists($paired_end{$name})){
   if($line[2]>${$paired_end{$name}}[1]){${$paired_end{$name}}[1]=$line[2];}
  }
  elsif(exists($paired_end{$line[4]}) && !(exists($paired_end{$name}))){
   ${$paired_end{$line[4]}}[0]=${$paired_end{$line[4]}}[0]+1;
   push(@{$paired_end{$line[4]}},$line[8]);
   if($line[2]>${$paired_end{$line[4]}}[2]){${$paired_end{$line[4]}}[2]=$line[2];}
   push(@{$paired_end{$name}},($line[1],$line[2],$line[7]));
  }
 }
 elsif($line[0] eq "F"){     ### left
  if(!(exists($single_left{$line[4]}))){push(@{$single_left{$line[4]}},(1,1,$line[2],"F",$line[8]));push(@{$single_left{$name}},($line[1],$line[2],$line[7]));}
  elsif(exists($single_left{$line[4]}) && exists($single_left{$name})){
   if($line[2]>${$single_left{$name}}[1]){${$single_left{$name}}[1]=$line[2];}
  }
  elsif(exists($single_left{$line[4]}) && !(exists($single_left{$name}))){
   ${$single_left{$line[4]}}[0]=${$single_left{$line[4]}}[0]+1;
   push(@{$single_left{$line[4]}},$line[8]);
   if($line[2]>${$single_left{$line[4]}}[2]){${$single_left{$line[4]}}[2]=$line[2];}
   push(@{$single_left{$name}},($line[1],$line[2],$line[7]));
  }
  if(exists($paired_end2{$line[4]}) && exists($single_right{$name})){
   if(exists($paired_end2{$name})){
    if(${$single_left{$name}}[1]+${$single_right{$name}}[1]>${$paired_end2{$name}}[1]){
     ${$paired_end2{$name}}[1]=${$single_left{$name}}[1]+${$single_right{$name}}[1];
    }
   }
   else{
    ${$paired_end2{$line[4]}}[0]=${$paired_end2{$line[4]}}[0]+1;
    push(@{$paired_end2{$line[4]}},$line[8]);
    push(@{$paired_end2{$name}},($line[1]."|".${$single_right{$name}}[0],${$single_left{$name}}[1]+${$single_right{$name}}[1],$line[7]."|".${$single_right{$name}}[2]));
    if(${$paired_end2{$name}}[1]>${$paired_end2{$line[4]}}[2]){${$paired_end2{$line[4]}}[2]=${$paired_end2{$name}}[1];}
   }
  }
  elsif(exists($single_right{$name}) && !(exists($paired_end2{$line[4]}))){
   push(@{$paired_end2{$name}},($line[1]."|".${$single_right{$name}}[0],${$single_left{$name}}[1]+${$single_right{$name}}[1],$line[7]."|".${$single_right{$name}}[2]));
   push(@{$paired_end2{$line[4]}},(1,1,${$single_left{$name}}[1]+${$single_right{$name}}[1],"P2",$line[8])); 
  }
 }
 elsif($line[0] eq "R"){     ### right
  if(!(exists($single_right{$line[4]}))){push(@{$single_right{$line[4]}},(1,1,$line[2],"R",$line[8]));push(@{$single_right{$name}},($line[1],$line[2],$line[7]));}
  elsif(exists($single_right{$line[4]}) && exists($single_right{$name})){
   if($line[2]>${$single_right{$name}}[1]){${$single_right{$name}}[1]=$line[2];}
  }
  elsif(exists($single_right{$line[4]}) && !(exists($single_right{$name}))){
   ${$single_right{$line[4]}}[0]=${$single_right{$line[4]}}[0]+1;
   push(@{$single_right{$line[4]}},$line[8]);
   if($line[2]>${$single_right{$line[4]}}[2]){${$single_right{$line[4]}}[2]=$line[2];}
   push(@{$single_right{$name}},($line[1],$line[2],$line[7]));
  }
  if(exists($paired_end2{$line[4]}) && exists($single_left{$name})){
   if(exists($paired_end2{$name})){ 
    if(${$single_left{$name}}[1]+${$single_right{$name}}[1]>${$paired_end2{$name}}[1]){
     ${$paired_end2{$name}}[1]=${$single_left{$name}}[1]+${$single_right{$name}}[1];  
    }
   }   
   else{
    ${$paired_end2{$line[4]}}[0]=${$paired_end2{$line[4]}}[0]+1;
    push(@{$paired_end2{$line[4]}},$line[8]);
    push(@{$paired_end2{$name}},(${$single_left{$name}}[0]."|".$line[1],${$single_left{$name}}[1]+${$single_right{$name}}[1],${$single_left{$name}}[2]."|".$line[7]));
    if(${$paired_end2{$name}}[1]>${$paired_end2{$line[4]}}[2]){${$paired_end2{$line[4]}}[2]=${$paired_end2{$name}}[1];}
   }   
  }
  elsif(exists($single_left{$name}) && !(exists($paired_end2{$line[4]}))){
   push(@{$paired_end2{$name}},(${$single_left{$name}}[0]."|".$line[1],${$single_left{$name}}[1]+${$single_right{$name}}[1],${$single_left{$name}}[2]."|".$line[7]));
   push(@{$paired_end2{$line[4]}},(1,1,${$single_left{$name}}[1]+${$single_right{$name}}[1],"P2",$line[8]));
  }
 }
                 }           ### 1. loop end

##################### for unique hash
my %unique=();
my %multiple=();
my %unique2=();
my @hash1=();
##################### mthod 1 2. loop start
my @keys1=keys %paired_end;
for(my $k=0;$k<@keys1;$k++){
 my @type=split /\|/, $keys1[$k];
 unless($#type eq 0){next;}
 $name=$keys1[$k]."|".${$paired_end{$keys1[$k]}}[4];
 if(!(exists($paired_end2{$keys1[$k]})) && !(exists($single_left{$keys1[$k]})) && !(exists($single_right{$keys1[$k]}))){
   if(${$paired_end{$keys1[$k]}}[0] eq 1){
    $unique2{$keys1[$k]}=0;
    @hash1=@{$paired_end{$keys1[$k]}};
    if(exists($unique{$hash1[4]})){${$unique{$hash1[4]}}[0]++;${$unique{$hash1[4]}}[1]++;}
    else{push(@{$unique{$hash1[4]}},(1,1,0,0,0,0,0));}           ##### total unique (P,P2,Single); Multiple (P,P2,Single);
    push(@{$unique{$hash1[4]}},($keys1[$k],@{$paired_end{$keys1[$k]}}[0..3],@{$paired_end{$name}}));
   }
  }
 }

################# paired_end2
my @keys1=keys %paired_end2;
for(my $k=0;$k<@keys1;$k++){
 my @type=split /\|/, $keys1[$k];
 unless($#type eq 0){next;}
 $name=$keys1[$k]."|".${$paired_end2{$keys1[$k]}}[4];
 if(!(exists($paired_end{$keys1[$k]}))){
   if(${$paired_end2{$keys1[$k]}}[0] eq 1){
    $unique2{$keys1[$k]}=0;
    @hash1=@{$paired_end2{$keys1[$k]}};
    if(exists($unique{$hash1[4]})){${$unique{$hash1[4]}}[0]++;${$unique{$hash1[4]}}[2]++;}
    else{push(@{$unique{$hash1[4]}},(1,0,1,0,0,0,0));}           ##### total unique (P,P2,Single); Multiple (P,P2,Single);
    push(@{$unique{$hash1[4]}},($keys1[$k],@{$paired_end2{$keys1[$k]}}[0..3],@{$paired_end2{$name}}));
   }
  }
 }

################ single_left
my @keys1=keys %single_left;
for(my $k=0;$k<@keys1;$k++){
 $name=$keys1[$k]."|".${$single_left{$keys1[$k]}}[4];
 my @type=split /\|/, $keys1[$k];
 unless($#type eq 0){next;}
 if(!(exists($paired_end{$keys1[$k]})) && !(exists($paired_end2{$keys1[$k]})) && !(exists($single_right{$line[4]}))){  
   if(${$single_left{$keys1[$k]}}[0] eq 1){
    $unique2{$keys1[$k]}=0;
    @hash1=@{$single_left{$keys1[$k]}};
    if(exists($unique{$hash1[4]})){${$unique{$hash1[4]}}[0]++;${$unique{$hash1[4]}}[3]++;}
    else{push(@{$unique{$hash1[4]}},(1,0,0,1,0,0,0));}           ##### total unique (P,P2,Single); Multiple (P,P2,Single);
    push(@{$unique{$hash1[4]}},($keys1[$k],@{$single_left{$keys1[$k]}}[0..3],@{$single_left{$name}}));
   }
  }
 }

################ single_right
my @keys1=keys %single_right;
for(my $k=0;$k<@keys1;$k++){
 my @type=split /\|/, $keys1[$k];
 unless($#type eq 0){next;}
 $name=$keys1[$k]."|".${$single_right{$keys1[$k]}}[4];
 if(!(exists($paired_end{$keys1[$k]})) && !(exists($paired_end2{$keys1[$k]})) && !(exists($single_left{$line[4]}))){
   if(${$single_right{$keys1[$k]}}[0] eq 1){
    $unique2{$keys1[$k]}=0;
    @hash1=@{$single_right{$keys1[$k]}};
    if(exists($unique{$hash1[4]})){${$unique{$hash1[4]}}[0]++;${$unique{$hash1[4]}}[3]++;}
    else{push(@{$unique{$hash1[4]}},(1,0,0,1,0,0,0));}           ##### total unique (P,P2,Single); Multiple (P,P2,Single);
    push(@{$unique{$hash1[4]}},($keys1[$k],@{$single_right{$keys1[$k]}}[0..3],@{$single_right{$name}}));
   }
  }
 }

my @temp1=keys %paired_end;
for(my $i=0;$i<@temp1;$i++){
 if(exists($unique2{$temp1[$i]})){next;}
 my @type=split /\|/, $temp1[$i];
 unless($#type eq 0){next;}
 my @temp2=@{$paired_end{$temp1[$i]}};
 my $candidate_virus="N";
 my $score=0;
 for(my $j=4;$j<@temp2;$j++){
  my $name2=$temp1[$i]."|".$temp2[$j];
  if(exists($unique{$temp2[$j]})){if(${$unique{$temp2[$j]}}[0]>$score){$score=${$unique{$temp2[$j]}}[0];$candidate_virus=$temp2[$j];}}
 }
 my $name2=$temp1[$i]."|".$candidate_virus;
 $unique2{$temp1[$i]}=0;
 if(exists($unique{$candidate_virus})){
  push(@{$unique{$candidate_virus}},($temp1[$i],@temp2[0..3],@{$paired_end{$name2}}));
  ${$unique{$candidate_virus}}[0]++;${$unique{$candidate_virus}}[4]++;
 }
}

########### for paired_end2
my @temp1=keys %paired_end2;
for(my $i=0;$i<@temp1;$i++){ 
 if(exists($unique2{$temp1[$i]})){next;}
 my @temp2=@{$paired_end2{$temp1[$i]}};
 my @type=split /\|/, $temp1[$i];
 unless($#type eq 0){next;}
 my $candidate_virus="N";
 my $score=0;
 for(my $j=4;$j<@temp2;$j++){
  my $name2=$temp1[$i]."|".$temp2[$j];
  if(exists($unique{$temp2[$j]})){if(${$unique{$temp2[$j]}}[0]>$score){$score=${$unique{$temp2[$j]}}[0];$candidate_virus=$temp2[$j];}}
 }
 $unique2{$temp1[$i]}=0;
 my $name2=$temp1[$i]."|".$candidate_virus;   
 if(exists($unique{$candidate_virus})){
  push(@{$unique{$candidate_virus}},($temp1[$i],@temp2[0..3],@{$paired_end2{$name2}}));
  ${$unique{$candidate_virus}}[0]++;${$unique{$candidate_virus}}[5]++;
 }
}

############ for single_left
my @temp1=keys %single_left;
for(my $i=0;$i<@temp1;$i++){
 if(exists($unique2{$temp1[$i]})){next;}
 my @temp2=@{$single_left{$temp1[$i]}};
 my @type=split /\|/, $temp1[$i];
 unless($#type eq 0){next;}
 my $candidate_virus="N";
 my $score=0;
 my @string=();
 for(my $j=4;$j<@temp2;$j++){
  my $name2=$temp1[$i]."|".$temp2[$j];
  if(exists($unique{$temp2[$j]})){if(${$unique{$temp2[$j]}}[0]>$score){@string=@{$single_left{$name2}};$score=${$unique{$temp2[$j]}}[0];$candidate_virus=$temp2[$j];}}
 }
 if(exists($single_right{$temp1[$i]})){
  @temp2=@{$single_right{$temp1[$i]}};
  for(my $j=4;$j<@temp2;$j++){
   my $name2=$temp1[$i]."|".$temp2[$j];
   if(exists($unique{$temp2[$j]})){if(${$unique{$temp2[$j]}}[0]>$score){@string=@{$single_right{$name2}};$score=${$unique{$temp2[$j]}}[0];$candidate_virus=$temp2[$j];}}
  }
 }
 my $name2=$temp1[$i]."|".$candidate_virus;
 $unique2{$temp1[$i]}=0;
 if(exists($unique{$candidate_virus})){
  push(@{$unique{$candidate_virus}},($temp1[$i],@temp2[0..3],@string));
  ${$unique{$candidate_virus}}[0]++;${$unique{$candidate_virus}}[6]++;
 }
}

############# for single_right
my @temp1=keys %single_right;
for(my $i=0;$i<@temp1;$i++){
 if(exists($unique2{$temp1[$i]})){next;}
 my @temp2=@{$single_right{$temp1[$i]}};
 my @type=split /\|/, $temp1[$i];
 unless($#type eq 0){next;}
 my $candidate_virus="N";
 my $score=0;
 my @string=();
 for(my $j=4;$j<@temp2;$j++){
  my $name2=$temp1[$i]."|".$temp2[$j];
  if(exists($unique{$temp2[$j]})){
   if(${$unique{$temp2[$j]}}[0]>$score){
     @string=@{$single_right{$name2}};
     $score=${$unique{$temp2[$j]}}[0];
     $candidate_virus=$temp2[$j];
   }
  }
 }
 if(exists($single_left{$temp1[$i]})){
  @temp2=@{$single_left{$temp1[$i]}};
  for(my $j=4;$j<@temp2;$j++){
   my $name2=$temp1[$i]."|".$temp2[$j];
   if(exists($unique{$temp2[$j]})){if(${$unique{$temp2[$j]}}[0]>$score){@string=@{$single_right{$name2}};$score=${$unique{$temp2[$j]}}[0];$candidate_virus=$temp2[$j];}}
  }
 } 
 $unique2{$temp1[$i]}=0;
 if(exists($unique{$candidate_virus})){
  push(@{$unique{$candidate_virus}},($temp1[$i],@temp2[0..3],@string));
  ${$unique{$candidate_virus}}[0]++;${$unique{$candidate_virus}}[6]++;
 }
}

my @temp1=keys %unique;
for(my $i=0;$i<@temp1;$i++){print"$temp1[$i] @{$unique{$temp1[$i]}}\n";}




