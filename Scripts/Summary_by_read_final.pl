#!usr/bin/perl
use strict;
my $line="";
my @line=();
my @reads=();
my $sample=$ARGV[0];
my %parameters=();
my %parameters2=();
my $name="";
my $i=0;

################################### open default_CS
open F1,"$sample.CS2";

while(<F1>){
 @reads=split;
 for($i=0;$i<@reads;$i+=159){
  @line=@reads[$i..$i+158];
  $name=join("\$",@line[0..15]);
  @{$parameters{$name}}=(0) x 12;
  # check virus: paired-end from target virus; target location or not; multiple map or not;
  # check human: paired end from human; target location or not; multiple map or not
  # check vicaller: paired-end from virus; target virus and location; multiple map on the target location; map to other viruses?
  # check homolog: between human and virus for each read

 if($line[8] eq "F" && $line[10] eq "S"){                ### F on human, R on virus
  if($i ne 477){
  if($line[26]>0) {${$parameters{$name}}[0]=-1;}           ### PE on virus
  elsif($line[29]>0) {${$parameters{$name}}[0]=-2;}
  elsif($line[32]>0) {${$parameters{$name}}[0]=-3;}
  else {${$parameters{$name}}[0]=1;}
  if($line[56]>0) {${$parameters{$name}}[1]=1;}                       ### target and multiple alignment on virus
  elsif($line[59]>0){${$parameters{$name}}[1]=2;}
  else{${$parameters{$name}}[1]=-1;}
  if($line[56]>0 && ($line[58] >= $line[55] || $line[61] >= $line[55] || $line[64] >= $line[55] || $line[67] >= $line[55])){${$parameters{$name}}[2]=-1;}
  elsif($line[59]>0 && ($line[55] >= $line[58] || $line[61] >= $line[58] || $line[64] >= $line[58] || $line[67] >= $line[58])){${$parameters{$name}}[2]=-2;}
  else{${$parameters{$name}}[2]=1;}

  if($line[116]>0) {${$parameters{$name}}[3]=-1;}
  elsif($line[119]>0) {${$parameters{$name}}[3]=-2;}
  elsif($line[122]>0) {${$parameters{$name}}[3]=-3;}
  else {${$parameters{$name}}[3]=1;}
  if($line[131]>0) {${$parameters{$name}}[4]=1;}
  elsif ($line[134]>0){${$parameters{$name}}[4]=2;}
  else{${$parameters{$name}}[4]=-1;}
  if($line[131]>0 && ($line[133] >= $line[130] || $line[136] >= $line[130] || $line[139] >= $line[130] || $line[142] >= $line[130])){${$parameters{$name}}[5]=-1;}
  elsif($line[134]>0 && ($line[130] >= $line[133] || $line[136] >= $line[133] || $line[139] >= $line[133] || $line[142] >= $line[133])){${$parameters{$name}}[5]=-2;}
  else{${$parameters{$name}}[5]=1;}

  if($line[71]>0) {${$parameters{$name}}[6]=-1;}
  elsif($line[74]>0) {${$parameters{$name}}[6]=-2;}
  elsif($line[77]>0) {${$parameters{$name}}[6]=-3;}
  else {${$parameters{$name}}[6]=1;}
  if($line[101]>0) {${$parameters{$name}}[7]=1;}
  elsif ($line[104]>0){${$parameters{$name}}[7]=2;}
  else{${$parameters{$name}}[7]=-1;}
  if($line[101]>0 && ($line[103] >= $line[100] || $line[106] >= $line[100] || $line[109] >= $line[100])){${$parameters{$name}}[8]=-1;}
  elsif($line[104]>0 && ($line[100] >= $line[103] || $line[106] >= $line[103] || $line[109] >= $line[103])){${$parameters{$name}}[8]=-2;}
  else{${$parameters{$name}}[8]=1;}
  if($line[101]>0 && $line[112] >= $line[100]){${$parameters{$name}}[9]=-1;}
  elsif($line[104]>0 && $line[112] >= $line[103]){${$parameters{$name}}[9]=-2;}
  else{${$parameters{$name}}[9]=1;}
  if($line[131]>0 && ($line[85] >= $line[130] || $line[88] >= $line[130] || $line[91] >= $line[130] || $line[94] >= $line[130] || $line[97] >= $line[130])){${$parameters{$name}}[10]=-1;}
  elsif($line[134]>0 && ($line[85] >= $line[133] || $line[88] >= $line[133] || $line[91] >= $line[133] || $line[94] >= $line[133] || $line[97] >= $line[133])){${$parameters{$name}}[10]=-2;}
  else{${$parameters{$name}}[10]=1;}
  if($line[56]>0 && ($line[145] >= $line[55] || $line[148] >= $line[55] || $line[151] >= $line[55] || $line[154] >= $line[55] || $line[157] >= $line[55])){${$parameters{$name}}[11]=-1;}
  elsif($line[59]>0 && ($line[145] >= $line[58] || $line[148] >= $line[58] || $line[151] >= $line[58] || $line[154] >= $line[58] || $line[157] >= $line[58])){${$parameters{$name}}[11]=-2;}
  else{${$parameters{$name}}[11]=1;}
 } else {                          ######## blast
  if($line[26]>0) {${$parameters{$name}}[0]=-1;}           ### PE on virus
  elsif($line[29]>0) {${$parameters{$name}}[0]=-2;}
  elsif($line[32]>0) {${$parameters{$name}}[0]=-3;}
  else {${$parameters{$name}}[0]=1;}
  if($line[56]>0) {${$parameters{$name}}[1]=1;}                       ### target and multiple alignment on virus
  elsif ($line[59]>0){${$parameters{$name}}[1]=2;}
  else{${$parameters{$name}}[1]=-1;}
  if($line[56]>0 && (($line[59] >0 && $line[58] <= $line[55]) || ($line[62] >0 && $line[61] <= $line[55]) || ($line[65] >0 && $line[64] <= $line[55]) || ($line[68] >0 && $line[67] <= $line[55]))){${$parameters{$name}}[2]=-1;}
   elsif($line[59]>0 && (($line[56] >0 && $line[55] <= $line[58]) || ($line[62] >0 && $line[61] <= $line[58]) || ($line[65] >0 && $line[64] <= $line[58]) || ($line[68] >0 && $line[67] <= $line[58]))){${$parameters{$name}}[2]=-2;}
   else{${$parameters{$name}}[2]=1;}
 
  if($line[116]>0) {${$parameters{$name}}[3]=-1;}
  elsif($line[119]>0) {${$parameters{$name}}[3]=-2;}
  elsif($line[122]>0) {${$parameters{$name}}[3]=-3;}
  else {${$parameters{$name}}[3]=1;}
  if($line[131]>0) {${$parameters{$name}}[4]=1;}
  elsif ($line[134]>0){${$parameters{$name}}[4]=2;}
  else{${$parameters{$name}}[4]=-1;}
  if($line[131]>0 && (($line[133] <= $line[130] && $line[134] > 0) || ($line[136] <= $line[130] && $line[137] > 0) || ($line[139] <= $line[130] && $line[140] > 0) || ($line[142] <= $line[130] && $line[143] > 0))){${$parameters{$name}}[5]=-1;}
  elsif($line[134]>0 && (($line[130] <= $line[133] && $line[131] > 0) || ($line[136] <= $line[133] && $line[137] > 0) || ($line[139] <= $line[133] && $line[140] > 0) || ($line[142] <= $line[133]&& $line[143] > 0))){${$parameters{$name}}[5]=-2;}
  else{${$parameters{$name}}[5]=1;}

  if($line[71]>0) {${$parameters{$name}}[6]=-1;}
  elsif($line[74]>0) {${$parameters{$name}}[6]=-2;}
  elsif($line[77]>0) {${$parameters{$name}}[6]=-3;}
  else {${$parameters{$name}}[6]=1;}
  if($line[101]>0) {${$parameters{$name}}[7]=1;}
  elsif ($line[104]>0){${$parameters{$name}}[7]=2;}
  else{${$parameters{$name}}[7]=-1;}
  if($line[101]>0 && (($line[104]>0 && $line[103] <= $line[100]) || ($line[107]>0 && $line[106] <= $line[100]) || ($line[110]>0 && $line[109] <= $line[100]))){${$parameters{$name}}[8]=-1;}
  elsif($line[104]>0 && (($line[101]>0 && $line[100] <= $line[103]) || ($line[107]>0 && $line[106] <= $line[103]) || ($line[110]>0 && $line[109] <= $line[103]))){${$parameters{$name}}[8]=-2;}
  else{${$parameters{$name}}[8]=1;}
  if($line[101]>0 && ($line[113]>0 && $line[112] <= $line[100])){${$parameters{$name}}[9]=-1;}
  elsif($line[104]>0 && ($line[113]>0 && $line[112] <= $line[103])){${$parameters{$name}}[9]=-2;}
  else{${$parameters{$name}}[9]=1;}
  if($line[131]>0 && (($line[86]>0 && $line[85] <= $line[130]) || ($line[89]>0 && $line[88] <= $line[130]) || ($line[92]>0 && $line[91] <= $line[130]) || ($line[95]>0 && $line[94] <= $line[130]) || ($line[98]>0 && $line[97] <= $line[130]))){${$parameters{$name}}[10]=-1;}
  elsif($line[134]>0 && (($line[86]>0 && $line[85] <= $line[133]) || ($line[89]>0 && $line[88] <= $line[133]) || ($line[92]>0 && $line[91] <= $line[133]) || ($line[95]>0 && $line[94] <= $line[133]) || ($line[98]>0 && $line[97] <= $line[133]))){${$parameters{$name}}[10]=-2;}
  else{${$parameters{$name}}[10]=1;}
  if($line[56]>0 && (($line[146]>0 && $line[145] <= $line[55]) || ($line[149]>0 && $line[148] <= $line[55]) || ($line[152]>0 && $line[151] <= $line[55]) || ($line[155]>0 && $line[154] <= $line[55]) || ($line[158]>0 && $line[157] <= $line[55]))){${$parameters{$name}}[11]=-1;}
  elsif($line[59]>0 && (($line[146]>0 && $line[145] <= $line[58]) || ($line[149]>0 && $line[148] <= $line[58]) || ($line[152]>0 && $line[151] <= $line[58]) || ($line[155]>0 && $line[154] <= $line[58]) || ($line[158]>0 && $line[157] <= $line[58]))){${$parameters{$name}}[11]=-2;}
  else{${$parameters{$name}}[11]=1;}
 }
}   ##### FR end  

#################### if second strand aligned to the human reference genome, and first strand aligned to the viral reference genome
 elsif($line[8] eq "S" && $line[10] eq "F"){
  if($i ne 477){             ##### bwa & blat
  if($line[26]>0) {${$parameters{$name}}[0]=-1;}           ### PE on virus
  elsif($line[29]>0) {${$parameters{$name}}[0]=-2;}
  elsif($line[32]>0) {${$parameters{$name}}[0]=-3;}
  else {${$parameters{$name}}[0]=1;}
  if($line[41]>0) {${$parameters{$name}}[1]=1;}                       ### target and multiple alignment on virus
  elsif ($line[44]>0){${$parameters{$name}}[1]=2;}
  else{${$parameters{$name}}[1]=-1;}
  if($line[41]>0 && ($line[43] >= $line[40] || $line[46] >= $line[40] || $line[49] >= $line[40] || $line[52] >= $line[40])){${$parameters{$name}}[2]=-1;}
  elsif($line[44]>0 && ($line[40] >= $line[43] || $line[46] >= $line[43] || $line[49] >= $line[43] || $line[52] >= $line[43])){${$parameters{$name}}[2]=-2;}
  else{${$parameters{$name}}[2]=1;}
 
  if($line[116]>0) {${$parameters{$name}}[3]=-1;}
  elsif($line[119]>0) {${$parameters{$name}}[3]=-2;}
  elsif($line[122]>0) {${$parameters{$name}}[3]=-3;}
  else {${$parameters{$name}}[3]=1;}
  if($line[146]>0) {${$parameters{$name}}[4]=1;}
  elsif ($line[149]>0){${$parameters{$name}}[4]=2;}
  else{${$parameters{$name}}[4]=-1;}
  if($line[146]>0 && ($line[148] >= $line[145] || $line[151] >= $line[145] || $line[154] >= $line[145] || $line[157] >= $line[145])){${$parameters{$name}}[5]=-1;}
  elsif($line[150]>0 && ($line[145] >= $line[148] || $line[151] >= $line[148] || $line[154] >= $line[148] || $line[157] >= $line[148])){${$parameters{$name}}[5]=-2;}
  else{${$parameters{$name}}[5]=1;}

  if($line[71]>0) {${$parameters{$name}}[6]=-1;}
  elsif($line[74]>0) {${$parameters{$name}}[6]=-2;}
  elsif($line[77]>0) {${$parameters{$name}}[6]=-3;}
  else {${$parameters{$name}}[6]=1;}
  if($line[86]>0) {${$parameters{$name}}[7]=1;}
  elsif ($line[89]>0){${$parameters{$name}}[7]=2;}
  else{${$parameters{$name}}[7]=-1;}
  if($line[86]>0 && ($line[88] >= $line[85] || $line[91] >= $line[85] || $line[94] >= $line[85])){${$parameters{$name}}[8]=-1;}
  elsif($line[89]>0 && ($line[85] >= $line[88] || $line[91] >= $line[88] || $line[94] >= $line[88])){${$parameters{$name}}[8]=-2;}
  else{${$parameters{$name}}[8]=1;}
  if($line[86]>0 && $line[97] >= $line[85]){${$parameters{$name}}[9]=-1;}
  elsif($line[89]>0 && $line[97] >= $line[88]){${$parameters{$name}}[9]=-2;}
  else{${$parameters{$name}}[9]=1;}

  if($line[146]>0 && ($line[100] >= $line[145] || $line[103] >= $line[145] || $line[106] >= $line[145] || $line[109] >= $line[145] || $line[112] >= $line[145])){${$parameters{$name}}[10]=-1;}
  elsif($line[149]>0 && ($line[100] >= $line[148] || $line[103] >= $line[148] || $line[106] >= $line[148] || $line[109] >= $line[148] || $line[112] >= $line[148])){${$parameters{$name}}[10]=-2;}
  else{${$parameters{$name}}[10]=1;}
  if($line[41]>0 && ($line[130] >= $line[40] || $line[133] >= $line[40] || $line[136] >= $line[40] || $line[139] >= $line[40] || $line[142] >= $line[40])){${$parameters{$name}}[11]=-1;}
  elsif($line[44]>0 && ($line[130] >= $line[43] || $line[133] >= $line[43] || $line[136] >= $line[43] || $line[139] >= $line[43] || $line[142] >= $line[43])){${$parameters{$name}}[11]=-2;}
  else{${$parameters{$name}}[11]=1;}
  }else {                           ##### blast
   if($line[26]>0) {${$parameters{$name}}[0]=-1;}           ### PE on virus
  elsif($line[29]>0) {${$parameters{$name}}[0]=-2;}
  elsif($line[32]>0) {${$parameters{$name}}[0]=-3;}
  else {${$parameters{$name}}[0]=1;}
  if($line[41]>0) {${$parameters{$name}}[1]=1;}                       ### target and multiple alignment on virus
  elsif ($line[44]>0){${$parameters{$name}}[1]=2;}
  else{${$parameters{$name}}[1]=-1;}
  if($line[41]>0 && (($line[44] > 0 && $line[43] <= $line[40]) || ($line[47] > 0 && $line[46] <= $line[40]) || ($line[50] > 0 && $line[49] <= $line[40]) || ($line[53] > 0 && $line[52] <= $line[40]))){${$parameters{$name}}[2]=-1;}
  elsif($line[44]>0 && (($line[41] > 0 && $line[40] <= $line[43]) || ($line[47] > 0 && $line[46] <= $line[43]) || ($line[50] > 0 && $line[49] <= $line[43]) || ($line[53] > 0 && $line[52] <= $line[43]))){${$parameters{$name}}[2]=-2;}
  else{${$parameters{$name}}[2]=1;}

  if($line[116]>0) {${$parameters{$name}}[3]=-1;}
  elsif($line[119]>0) {${$parameters{$name}}[3]=-2;}
  elsif($line[122]>0) {${$parameters{$name}}[3]=-3;}
  else {${$parameters{$name}}[3]=1;}
  if($line[146]>0) {${$parameters{$name}}[4]=1;}
  elsif ($line[149]>0){${$parameters{$name}}[4]=2;}
  else{${$parameters{$name}}[4]=-1;}
  if($line[146]>0 && (($line[149] > 0 && $line[148] <= $line[145]) || ($line[152] > 0 && $line[151] <= $line[145]) || ($line[155] > 0 && $line[154] <= $line[145]) || ($line[158] > 0 && $line[157] <= $line[145]))){${$parameters{$name}}[5]=-1;}
  elsif($line[149]>0 && (($line[146] > 0 && $line[145] <= $line[148]) || ($line[152] > 0 && $line[151] <= $line[148]) || ($line[155] > 0 && $line[154] <= $line[148]) || ($line[158] > 0 && $line[157] <= $line[148]))){${$parameters{$name}}[5]=-2;}
  else{${$parameters{$name}}[5]=1;}

  if($line[71]>0) {${$parameters{$name}}[6]=-1;}
  elsif($line[74]>0) {${$parameters{$name}}[6]=-2;}
  elsif($line[77]>0) {${$parameters{$name}}[6]=-3;}
  else {${$parameters{$name}}[6]=1;}
  if($line[86]>0) {${$parameters{$name}}[7]=1;}
  elsif ($line[91]>0){${$parameters{$name}}[7]=2;}
  else{${$parameters{$name}}[7]=-1;}
  if($line[86]>0 && (($line[89] > 0 && $line[88] <= $line[85]) || ($line[92] > 0 && $line[91] <= $line[85]) || ($line[95] > 0 && $line[94] <= $line[85]))){${$parameters{$name}}[8]=-1;}
  elsif($line[89]>0 && (($line[86] > 0 && $line[85] <= $line[88]) || ($line[92] > 0 && $line[91] <= $line[88]) || ($line[95] > 0 && $line[94] <= $line[88]))){${$parameters{$name}}[8]=-2;}
  else{${$parameters{$name}}[8]=1;}
  if($line[86]>0 && $line[98] > 0 && $line[97] <= $line[85]){${$parameters{$name}}[9]=-1;}
  elsif($line[89]>0 && $line[98] > 0 && $line[97] <= $line[88]){${$parameters{$name}}[9]=-2;}
  else{${$parameters{$name}}[9]=1;}

  if($line[146]>0 && (($line[101] > 0 && $line[100] <= $line[145]) || ($line[104] > 0 && $line[103] <= $line[145]) || ($line[107] > 0 && $line[106] <= $line[145]) || ($line[110] > 0 && $line[109] <= $line[145]) || ($line[113] > 0 && $line[112] <= $line[145]))){${$parameters{$name}}[10]=-1;}
  elsif($line[149]>0 && (($line[101] > 0 && $line[100] <= $line[148]) || ($line[104] > 0 && $line[103] <= $line[148]) || ($line[107] > 0 && $line[106] <= $line[148]) || ($line[110] > 0 && $line[109] <= $line[148]) || ($line[113] > 0 && $line[112] <= $line[148]))){${$parameters{$name}}[10]=-2;}
  else{${$parameters{$name}}[10]=1;}
  if($line[41]>0 && (($line[131] > 0 && $line[130] <= $line[40]) || ($line[134] > 0 && $line[133] <= $line[40]) || ($line[137] > 0 && $line[136] <= $line[40]) || ($line[140] > 0 && $line[139] <= $line[40]) || ($line[143] > 0 && $line[142] <= $line[40]))){${$parameters{$name}}[11]=-1;}
  elsif($line[44]>0 && (($line[131] > 0 && $line[130] <= $line[43]) || ($line[134] > 0 && $line[133] <= $line[43]) || ($line[137] > 0 && $line[136] <= $line[43]) || ($line[140] > 0 && $line[139] <= $line[43]) || ($line[143] > 0 && $line[142] <= $line[43]))){${$parameters{$name}}[11]=-2;}
  else{${$parameters{$name}}[11]=1;}
}
 }    #### RF end

#################### if it is split reads
 elsif($line[8] eq "S1") {
  if($i ne 477){
  if($line[26]>0) {${$parameters{$name}}[0]=-1;}           ### PE on virus
  elsif($line[29]>0) {${$parameters{$name}}[0]=-2;}
  elsif($line[32]>0) {${$parameters{$name}}[0]=-3;}
  else {${$parameters{$name}}[0]=1;}
  if($line[41]>0) {${$parameters{$name}}[1]=1;}                       ### target and multiple alignment on virus
  elsif ($line[44]>0){${$parameters{$name}}[1]=2;}
  else{${$parameters{$name}}[1]=-1;}
  if($line[41]>0 && ($line[43] >= $line[40] || $line[46] >= $line[40] || $line[49] >= $line[40] || $line[52] >= $line[40])){${$parameters{$name}}[2]=-1;}
  elsif($line[44]>0 && ($line[40] >= $line[43] || $line[46] >= $line[43] || $line[49] >= $line[43] || $line[52] >= $line[43])){${$parameters{$name}}[2]=-2;}
  else{${$parameters{$name}}[2]=1;}

  if($line[71]>0) {${$parameters{$name}}[6]=-1;}
  elsif($line[74]>0) {${$parameters{$name}}[6]=-2;}
  elsif($line[77]>0) {${$parameters{$name}}[6]=-3;}
  else {${$parameters{$name}}[6]=1;}
  if($line[86]>0) {${$parameters{$name}}[7]=1;}
  elsif ($line[89]>0){${$parameters{$name}}[7]=2;}
  else{${$parameters{$name}}[7]=-1;}
  if($line[86]>0 && ($line[88] >= $line[85] || $line[91] >= $line[85] || $line[94] >= $line[85])){${$parameters{$name}}[8]=-1;}
  elsif($line[89]>0 && ($line[85] >= $line[88] || $line[91] >= $line[88] || $line[94] >= $line[88])){${$parameters{$name}}[8]=-2;}
  else{${$parameters{$name}}[8]=1;}
  if($line[86]>0 && $line[97] >= $line[85]){${$parameters{$name}}[9]=-1;}
  elsif($line[89]>0 && $line[97] >= $line[88]){${$parameters{$name}}[9]=-2;}
  else{${$parameters{$name}}[9]=1;}

  if($line[41]>0 && ($line[130] >= $line[40] || $line[133] >= $line[40] || $line[136] >= $line[40] || $line[139] >= $line[40] || $line[142] >= $line[40])){${$parameters{$name}}[11]=-1;}
  elsif($line[44]>0 && ($line[130] >= $line[43] || $line[133] >= $line[43] || $line[136] >= $line[43] || $line[139] >= $line[43] || $line[142] >= $line[43])){${$parameters{$name}}[11]=-2;}
  else{${$parameters{$name}}[11]=1;}
  } else {                    ##### blast
  if($line[26]>0) {${$parameters{$name}}[0]=-1;}           ### PE on virus  
  elsif($line[29]>0) {${$parameters{$name}}[0]=-2;}
  elsif($line[32]>0) {${$parameters{$name}}[0]=-3;}
  else {${$parameters{$name}}[0]=1;}
  if($line[41]>0) {${$parameters{$name}}[1]=1;}                       ### target and multiple alignment on virus  elsif ($line[45]>0){${$parameters{$name}}[1]=2;}  else{${$parameters{$name}}[1]=-1;}
  if($line[41]>0 && (($line[44] > 0 && $line[43] <= $line[40]) || ($line[47] > 0 && $line[46] <= $line[40]) || ($line[50] > 0 && $line[49] <= $line[40]) || ($line[53] > 0 && $line[52] <= $line[40]))){${$parameters{$name}}[2]=-1;}
  elsif($line[44]>0 && (($line[41] > 0 && $line[40] <= $line[43]) || ($line[47] > 0 && $line[46] <= $line[43]) || ($line[50] > 0 && $line[49] <= $line[43]) || ($line[53] > 0 && $line[52] <= $line[43]))){${$parameters{$name}}[2]=-2;}
  else{${$parameters{$name}}[2]=1;}
  if($line[71]>0) {${$parameters{$name}}[6]=-1;}
  elsif($line[74]>0) {${$parameters{$name}}[6]=-2;}
  elsif($line[77]>0) {${$parameters{$name}}[6]=-3;}
  else {${$parameters{$name}}[6]=1;}
  if($line[86]>0) {${$parameters{$name}}[7]=1;}
  elsif ($line[89]>0){${$parameters{$name}}[7]=2;}
  else{${$parameters{$name}}[7]=-1;}
  if($line[86]>0 && (($line[89] > 0 && $line[88] <= $line[85]) || ($line[92] > 0 && $line[91] <= $line[85]) || ($line[95] > 0 && $line[94] <= $line[85]))){${$parameters{$name}}[8]=-1;}
  elsif($line[89]>0 && (($line[86] > 0 && $line[85] <= $line[88]) || ($line[92] > 0 && $line[91] <= $line[88]) || ($line[95] > 0 && $line[94] <= $line[88]))){${$parameters{$name}}[8]=-2;}
  else{${$parameters{$name}}[8]=1;}
  if($line[86]>0 && $line[98] > 0 && $line[97] <= $line[85]){${$parameters{$name}}[9]=-1;}
  elsif($line[89]>0 && $line[98] > 0 && $line[97] <= $line[88]){${$parameters{$name}}[9]=-2;}
  else{${$parameters{$name}}[9]=1;}
  if($line[41]>0 && (($line[131] > 0 && $line[130] <= $line[40]) || ($line[134] > 0 && $line[133] <= $line[40]) || ($line[137] > 0 && $line[136] <= $line[40]) || ($line[140] > 0 && $line[139] <= $line[40]) || ($line[143] > 0 && $line[142] <= $line[40]))){${$parameters{$name}}[11]=-1;}
  elsif($line[44]>0 && (($line[131] > 0 && $line[130] <= $line[43]) || ($line[134] > 0 && $line[133] <= $line[43]) || ($line[137] > 0 && $line[136] <= $line[43]) || ($line[140] > 0 && $line[139] <= $line[43]) || ($line[143] > 0 && $line[142] <= $line[43]))){${$parameters{$name}}[11]=-2;}
  else{${$parameters{$name}}[11]=1;}
  }
 }   ##### S1 end
  push(@{$parameters2{$name}},@{$parameters{$name}});
  @{$parameters{$name}}=(0) x12;
 }
  push(@{$parameters2{$name}},@reads);
}               ##### while end

####################### Summary information of each read
my @read_name=keys %parameters2;
my $chimeric_number=0;
my $chimeric_correct=0;
my $chimeric_weak=0;
my $chimeric_wrong=0;
my $split_number=0;
my $split_correct=0;
my $split_weak=0;
my $split_wrong=0;

for(my $i=0;$i<@read_name;$i++){
 my @line2=@{$parameters2{$read_name[$i]}};
 if($line2[56] eq "S1"){
  $split_number++;
  if($line2[1]>=0 && $line2[7]>=0 && $line2[11]>=0  &&  $line2[13]>=0 && $line2[19]>=0 && $line2[23]>=0  &&  $line2[25]>=0 && $line2[31]>=0 && $line2[35]>=0  &&  $line2[37]>=0 && $line2[43]>=0 && $line2[47]>=0){$split_correct++;}
  elsif($line2[11]<0){$split_wrong++;}
  elsif($line2[11]<0 && $line2[23]<0 && $line2[35]<0 && $line2[47]<0){$split_wrong++;}
  elsif($line2[1]<0 && $line2[13]<0 && $line2[25]<0 && $line2[37]<0){$split_wrong++;}
  elsif($line2[7]<0 && $line2[19]<0 && $line2[31]<0 && $line2[43]<0){$split_wrong++;}
  else{$split_weak++;}
                       }
 else{
  $chimeric_number++; 
  if($line2[0]>=0 && $line2[1]>=0 &&$line2[3]>=0 &&$line2[4]>=0 &&$line2[5]>=0 && $line2[6]>=0 &&$line2[7]>=0 &&$line2[10]>=0 &&$line2[11]>=0 && $line2[12]>=0 &&$line2[13]>=0 &&$line2[15]>=0 &&$line2[16]>=0 &&$line2[17]>=0 &&$line2[18]>=0 &&$line2[19]>=0 &&$line2[22]>=0 &&$line2[23]>=0 &&$line2[24]>=0 &&$line2[25]>=0 &&$line2[27]>=0 &&$line2[28]>=0 &&$line2[29]>=0 &&$line2[30]>=0 &&$line2[31]>=0 && $line2[34]>=0 &&$line2[35]>=0 &&$line2[36]>=0 &&$line2[37]>=0 &&$line2[39]>=0 &&$line2[40]>=0&&$line2[41]>=0&&$line2[42]>=0&&$line2[43]>=0&&$line2[46]>=0 && $line2[47]>=0) {$chimeric_correct++;}
  elsif($line2[1]<0 && $line2[13]<0 && $line2[25]<0 && $line2[37]<0){$chimeric_wrong++;} ##### if cannot aligned to the target
  elsif($line2[24]<0 && $line2[36]<0){$chimeric_wrong++;}          			 ##### if both end map to virus
  elsif($line2[10]<0 && $line2[22]<0){$chimeric_wrong++;}				 ##### if one end map to both human and viral  genome
  elsif($line2[27]<0 && $line2[39]<0){$chimeric_wrong++;} 			         ##### if both end map to human
  else{$chimeric_weak++;}
 }
}


for(my $i=0;$i<@read_name;$i++){
 my @line2=@{$parameters2{$read_name[$i]}};
 print "$#read_name $chimeric_number $chimeric_correct $chimeric_weak $chimeric_wrong $split_number $split_correct $split_weak $split_wrong @line2\n";
}


