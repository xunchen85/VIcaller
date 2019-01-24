#!usr/bin/perl
use strict;

my $line="";
my @line=();
my %read=();
my $sample=$ARGV[0];

################################### open default_CS
open F1,"$sample.default_CS";
while(<F1>){
 @line=split;
 my $name=join(" ",@line[0..13]);
 @{$read{$name}}=@line;
}

################################## open notdefault_CS
open F1,"$sample.notdefault_CS";
while(<F1>){
 @line=split;
 my $name=join(" ",@line[0..13]);
 push(@{$read{$name}},@line);
}

################################## open psl_CS
open F1,"$sample.psl_CS";
while(<F1>){
 @line=split;
  my $name=join(" ",@line[0..13]);
   push(@{$read{$name}},@line);
   }

################################## open out_CS
open F1,"$sample.out_CS";
while(<F1>){
 @line=split;
  my $name=join(" ",@line[0..13]);
   push(@{$read{$name}},@line);
   }

##################################
my @temp1=keys %read;
for(my $i=0;$i<@temp1;$i++){
}

########################################################################
my $col=0;
my @results=(0) x 30;
my @read=();
for(my $i=0;$i<@temp1;$i++){
 $col=0;
 @read=@{$read{$temp1[$i]}};
 for(my $j=0;$j<@read;$j+=159){

##########
  if($read[$j+25] eq 0 && $read[$j+28] eq 0 && $read[$j+31] eq 0 && $read[$j+34] eq 0 && $read[$j+37] eq 0){$results[$col++]=0;}
  else{$results[$col++]=1;}
  if($read[$j+40] > 0 || $read[$j+43] > 0){$results[$col++]=1;}
  else{$results[$col++]=0;}
  if($read[$j+55] > 0 || $read[$j+58] > 0){$results[$col++]=1;}
  else{$results[$col++]=0;}

##########
  if($read[$j+70] eq 0 && $read[$j+73] eq 0 && $read[$j+76] eq 0 && $read[$j+79] eq 0 && $read[$j+82] eq 0){$results[$col++]=0;}
  else{$results[$col++]=1;$results[$col++]=1;}
  if($read[$j+85] > 0 || $read[$j+88] > 0){$results[$col++]=1;}
  else{$results[$col++]=0;}
  if($read[$j+100] > 0 || $read[$j+103] > 0){$results[$col++]=1;}
  else{$results[$col++]=0;}

##########
  if($read[$j+115] eq 0 && $read[$j+118] eq 0 && $read[$j+121] eq 0 && $read[$j+124] eq 0 && $read[$j+127] eq 0){$results[$col++]=0;}
  else{$results[$col++]=1;}
  if($read[$j+130] > 0 || $read[$j+133] > 0){$results[$col++]=1;}
  else{$results[$col++]=0;}
  if($read[$j+145] > 0 || $read[$j+148] > 0){$results[$col++]=1;} 
  else{$results[$col++]=0;}
 }

##########
my @results2=("-") x 4;
my $j=0;

for(my $i=0;$i<@results;$i+=9){
 if($read[8] eq "S" && $read[10] eq "F"){
  if($results[$i+0] eq 0 && $results[$i+1] eq 1 && $results[$i+2] eq 0){$results2[$j]=$results2[$j]."T";}        ##### validated in the target database
  if($results[$i+3] eq 0 && $results[$i+4] eq 1 && $results[$i+5] eq 0){$results2[$j]=$results2[$j]."V";}        ##### validated in the vicaller database
  if($results[$i+6] eq 0 && $results[$i+7] eq 0 && $results[$i+8] eq 1){$results2[$j]=$results2[$j]."H";}        ##### validated in the human reference genome
  }
 elsif($read[8] eq "F" && $read[10] eq "S"){
  if($results[$i+0] eq 0 && $results[$i+1] eq 0 && $results[$i+2] eq 1){$results2[$j]=$results2[$j]."T";}
  if($results[$i+3] eq 0 && $results[$i+4] eq 0 && $results[$i+5] eq 1){$results2[$j]=$results2[$j]."V";}
  if($results[$i+6] eq 0 && $results[$i+7] eq 1 && $results[$i+8] eq 0){$results2[$j]=$results2[$j]."H";}
  }
 elsif($read[8] eq "S1"){
  if($results[$i+0] eq 0 && $results[$i+1] eq 1 && $results[$i+2] eq 0){$results2[$j]=$results2[$j]."T";}
  if($results[$i+3] eq 0 && $results[$i+4] eq 1 && $results[$i+5] eq 0){$results2[$j]=$results2[$j]."V";}
  if($results[$i+6] eq 0 && $results[$i+7] eq 0 && $results[$i+8] eq 0){$results2[$j]=$results2[$j]."H";}
  }
 $j++;
}
print "@read[0..23] @results2\n";

##########
 @results=();
 @results=(0) x 30;
 @results2=();
 @results2=(0) x 4;
}

