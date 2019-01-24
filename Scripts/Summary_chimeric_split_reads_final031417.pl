#!usr/bin/perl
use strict;

############## variables, and names of input files;
my $line="";
my @line=();
my %read_id=();
my $name="";
my $type="";
my $gi="";                 #### define the target GI;

############## extract read seq
my %file=();
my $b1="";
my $read_name="";
my $temp04="";
my $temp05="";
my $file_name="";
my $aa="";
my @temp=();
my $read="";

open F1,"$ARGV[0]_aligned_1.fuq"||di("di");
while($read_name=<F1>){
   $read=$read_name;
   my @temp11=split /\s+/, $read;
   $read_name=$temp11[0];
   $read_name=~s/\/[12]\s+$//;
   $read_name=~s/^@//;
   $read_name=~s/\/1//;
   $read_name=~s/\/2//;
   my $temp01=<F1>;
   chmod($temp01);
   $temp01=~s/\s+$//;
   my $temp02=<F1>;
   my $temp03=<F1>;
   $temp04=$temp01;
   $file{$read_name}=$temp04;
#   print "$read_name\n";
}
close F1;

my %file2=();
open F1,"$ARGV[0]_aligned_2.fuq"||di("di");
while($read_name=<F1>){
   $read=$read_name;
   my @temp11=split /\s+/, $read;
   $read_name=$temp11[0];
   $read_name=~s/\/[12]\s+$//;
   $read_name=~s/^@//;
   $read_name=~s/\/1//;
   $read_name=~s/\/2//;
   my $temp01=<F1>;
   chmod($temp01);
   $temp01=~s/\s+$//;
   my $temp02=<F1>;
   chmod($temp02);
   my $temp03=<F1>;
   chmod($temp03);
   $temp04=$temp01;
   $file2{$read_name}=$temp04;
}
close F1;

open F1,"$ARGV[0]_aligned_sf.fuq"||di("di");
while($read_name=<F1>){
   $read=$read_name;
   my @temp11=split /\s+/, $read;
   $read_name=$temp11[0];
   $read_name=~s/\/[12]\s+$//;
   $read_name=~s/^@//;
   $read_name=~s/\/1//;
   $read_name=~s/\/2//;
   my $temp01=<F1>;
   chmod($temp01);
   $temp01=~s/\s+$//;
   my $temp02=<F1>;
   my $temp03=<F1>;
   $temp04=$temp01;
   $file{$read_name}=$temp04;
}
close F1;

################################# read the file containing full information of each chimeric and split read;
my @infor1=();
my @infor2=();
my %infor=();
open F1,"$ARGV[0].information";
while(<F1>){
 @infor1=split;
 ##################### 1, check if 0,1, Target; 2,3, Nearby region; 4,5, Other position; 6,7, Other GI; 8,9, Other Virus; 
 push(@infor1,(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0));          #### P
 push(@infor1,(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0));          #### F
 push(@infor1,(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0));          #### R
 push(@{$infor{$infor1[7]}},[@infor1]);
}

############## check which tool it had been used, and specified the column of start position of the read;
my $posi=0;
if($ARGV[1] =~ "_infection"){$posi=7;}
elsif($ARGV[1] =~ "\.out"){$posi=13;}                           ##### maybe will be 1 offset, because of different tools;
elsif($ARGV[1] =~ "\.psl"){$posi=20;}                           ##### maybe will be 1 offset, because of different tools;

############## check which database it had been aligned to;
my $db="";
my $col=0;

############## VIcaller
if($ARGV[1] =~ "vicaller"){  ##### start 1 of vicaller database
 $db="vicaller";
 open F2,"$ARGV[1]\n";
 while(<F2>){
  @line=split;
  @infor2=();
  my @read_type=split /\|/, $line[4];
  if($#read_type>1){$line[1]="F";}

######## which kind of read, P, PF/F, or PR/R;
  if($line[1] eq "P"){$col=14;}
  elsif($line[1] =~ "F"){$col=29;}
  elsif($line[1] =~ "R"){$col=44;}
  if(exists($infor{$line[4]})){
   for(my $i=0;$i<@{$infor{$line[4]}};$i++){
    @infor2=@{$infor{$line[4]}[$i]};
   if($infor2[4] eq $line[0] && $infor2[5] eq $line[2] && abs($line[$posi]-$infor2[13])<=1 && $line[3] >= $infor2[$col+1]){
    $infor2[$col+1]=$line[3];
    $infor2[$col]=$line[$posi];
    if($infor2[$col+2] eq 0){$infor2[$col+2]=1;}
    elsif($line[3] eq $infor2[$col+1]){$infor2[$col+2]++;}
    elsif($line[3] > $infor2[$col+1]){$infor2[$col+2]=1;}
   }
   elsif($infor2[4] eq $line[0] && $infor2[5] eq $line[2] && abs($line[$posi]-$infor2[13])>1 && abs($line[$posi]-$infor2[13])<=1000 && $line[3] >= $infor2[$col+4]){
    $infor2[$col+4]=$line[3];
    $infor2[$col+3]=$line[$posi];
    if($infor2[$col+5] eq 0){$infor2[$col+5]=1;}
    elsif($line[3] eq $infor2[$col+4]){$infor2[$col+5]++;}
    elsif($line[3] > $infor2[$col+4]){$infor2[$col+5]=1;}

   }
   elsif($infor2[4] eq $line[0] && $infor2[5] eq $line[2] && abs($line[$posi]-$infor2[13])>1000 && $line[3] >= $infor2[$col+7]){
    $infor2[$col+7]=$line[3];
    $infor2[$col+6]=$line[$posi];
    if($infor2[$col+8] eq 0){$infor2[$col+8]=1;}
    elsif($line[3] eq $infor2[$col+7]){$infor2[$col+8]++;}
    elsif($line[3] > $infor2[$col+7]){$infor2[$col+8]=1;}

   }
   elsif($infor2[4] eq $line[0] && $infor2[5] ne $line[2] && $line[3] >= $infor2[$col+10]){
    $infor2[$col+10]=$line[3];
    $infor2[$col+9]=$line[0]."|".$line[2]."|".$line[$posi];
    if($infor2[$col+11] eq 0){$infor2[$col+11]=1;}
    elsif($line[3] eq $infor2[$col+10]){$infor2[$col+11]++;}
    elsif($line[3] > $infor2[$col+10]){$infor2[$col+11]=1;}

   }
   elsif($infor2[4] ne $line[0] && $line[3] >= $infor2[$col+13]){
    $infor2[$col+13]=$line[3];
    $infor2[$col+12]=$line[0]."|".$line[2]."|".$line[$posi];
    if($infor2[$col+14] eq 0){$infor2[$col+14]=1;}
    elsif($line[3] eq $infor2[$col+13]){$infor2[$col+14]++;}
    elsif($line[3] > $infor2[$col+13]){$infor2[$col+14]=1;}

   }
   @{$infor{$line[4]}[$i]}=@infor2;
                                         }
                            }
           }
 close F2;
}                            ##### end 1

######################## target
elsif($ARGV[1] =~ "target"){ ##### start 2 of target reference genome
 $db="target";               
 open F2,"$ARGV[1]\n";
 while(<F2>){
  @line=split;
  my @read_type=split /\|/, $line[4];
  if($#read_type>1){$line[1]="F";}
  @infor2=();
  if($line[1] eq "P"){$col=14;}
  elsif($line[1] =~ "F"){$col=29;}
  elsif($line[1] =~ "R"){$col=44;}
   if(exists($infor{$line[4]})){
   for(my $i=0;$i<@{$infor{$line[4]}};$i++){
    @infor2=@{$infor{$line[4]}[$i]};
   if($infor2[4] eq $line[0] && $infor2[5] eq $line[2] && abs($line[$posi]-$infor2[13])<=1 && $line[3] >= $infor2[$col+1]){
    $infor2[$col+1]=$line[3];
    $infor2[$col]=$line[$posi];
    if($infor2[$col+2] eq 0){$infor2[$col+2]=1;}
    elsif($line[3] eq $infor2[$col+1]){$infor2[$col+2]++;}
    elsif($line[3] > $infor2[$col+1]){$infor2[$col+2]=1;}
   }
   elsif($infor2[4] eq $line[0] && $infor2[5] eq $line[2] && abs($line[$posi]-$infor2[13])>1 && abs($line[$posi]-$infor2[13])<=1000 && $line[3] >= $infor2[$col+4]){
    $infor2[$col+4]=$line[3];
    $infor2[$col+3]=$line[$posi];
    if($infor2[$col+5] eq 0){$infor2[$col+5]=1;}
    elsif($line[3] eq $infor2[$col+4]){$infor2[$col+5]++;}
    elsif($line[3] > $infor2[$col+4]){$infor2[$col+5]=1;}

   }
   elsif($infor2[4] eq $line[0] && $infor2[5] eq $line[2] && abs($line[$posi]-$infor2[13])>1000 && $line[3] >= $infor2[$col+7]){
    $infor2[$col+7]=$line[3];
    $infor2[$col+6]=$line[$posi];
    if($infor2[$col+8] eq 0){$infor2[$col+8]=1;}
    elsif($line[3] eq $infor2[$col+7]){$infor2[$col+8]++;}
    elsif($line[3] > $infor2[$col+7]){$infor2[$col+8]=1;}

   }
   elsif($infor2[4] eq $line[0] && $infor2[5] ne $line[2] && $line[3] >= $infor2[$col+10]){
    $infor2[$col+10]=$line[3];
    $infor2[$col+9]=$line[0]."|".$line[2]."|".$line[$posi];
    if($infor2[$col+11] eq 0){$infor2[$col+11]=1;}
    elsif($line[3] eq $infor2[$col+10]){$infor2[$col+11]++;}
    elsif($line[3] > $infor2[$col+10]){$infor2[$col+11]=1;}

   }
   elsif($infor2[4] ne $line[0] && $line[3] >= $infor2[$col+13]){
    $infor2[$col+13]=$line[3];
    $infor2[$col+12]=$line[0]."|".$line[2]."|".$line[$posi];
    if($infor2[$col+14] eq 0){$infor2[$col+14]=1;}
    elsif($line[3] eq $infor2[$col+13]){$infor2[$col+14]++;}
    elsif($line[3] > $infor2[$col+13]){$infor2[$col+14]=1;}

   }
   @{$infor{$line[4]}[$i]}=@infor2;
                                         }
                            }
           }
 close F2;
}                         ##### end 

######################## human
elsif($ARGV[1] =~ "human"){  ##### start3 of human reference genome
 $db="human";
 open F2,"$ARGV[1]\n";
 while(<F2>){
  @line=split;
  @infor2=();
  my @read_type=split /\|/, $line[4];
  if($#read_type>1){$line[1]="F";}
  if($line[1] eq "P"){$col=14;}
  elsif($line[1] =~ "F"){$col=29;}
  elsif($line[1] =~ "R"){$col=44;}
  if(exists($infor{$line[4]})){
   for(my $i=0;$i<@{$infor{$line[4]}};$i++){
    @infor2=@{$infor{$line[4]}[$i]};
   if($infor2[1] eq $line[0] && $infor2[1] eq $line[2] && abs($line[$posi]-$infor2[12])<=1 && $line[3] >= $infor2[$col+1]){
    $infor2[$col+1]=$line[3];
    $infor2[$col]=$line[$posi];
    if($infor2[$col+2] eq 0){$infor2[$col+2]=1;}
    elsif($line[3] eq $infor2[$col+1]){$infor2[$col+2]++;}
    elsif($line[3] > $infor2[$col+1]){$infor2[$col+2]=1;}

   }
   elsif($infor2[1] eq $line[0] && $infor2[1] eq $line[2] && abs($line[$posi]-$infor2[12])>1 && abs($line[$posi]-$infor2[12])<=1000 && $line[3] >= $infor2[$col+4]){
    $infor2[$col+4]=$line[3];
    $infor2[$col+3]=$line[$posi];
    if($infor2[$col+5] eq 0){$infor2[$col+5]=1;}
    elsif($line[3] eq $infor2[$col+4]){$infor2[$col+5]++;}
    elsif($line[3] > $infor2[$col+4]){$infor2[$col+5]=1;}

   }
   elsif($infor2[1] eq $line[0] && $infor2[1] eq $line[2] && abs($line[$posi]-$infor2[12])>1000 && $line[3] >= $infor2[$col+7]){
    $infor2[$col+7]=$line[3];
    $infor2[$col+6]=$line[$posi];
    if($infor2[$col+8] eq 0){$infor2[$col+8]=1;}
    elsif($line[3] eq $infor2[$col+7]){$infor2[$col+8]++;}
    elsif($line[3] > $infor2[$col+7]){$infor2[$col+8]=1;}

   }
   elsif($infor2[1] eq $line[0] && $infor2[1] ne $line[2] && $line[3] >= $infor2[$col+10]){
    $infor2[$col+10]=$line[3];
    $infor2[$col+9]=$line[0]."|".$line[2]."|".$line[$posi];
    if($infor2[$col+11] eq 0){$infor2[$col+11]=1;}
    elsif($line[3] eq $infor2[$col+10]){$infor2[$col+11]++;}
    elsif($line[3] > $infor2[$col+10]){$infor2[$col+11]=1;}

   }
   elsif($infor2[1] ne $line[0] && $line[3] >= $infor2[$col+13]){
    $infor2[$col+13]=$line[3];
    $infor2[$col+12]=$line[0]."|".$line[2]."|".$line[$posi];
    if($infor2[$col+14] eq 0){$infor2[$col+14]=1;}
    elsif($line[3] eq $infor2[$col+13]){$infor2[$col+14]++;}
    elsif($line[3] > $infor2[$col+13]){$infor2[$col+14]=1;}
   }
   @{$infor{$line[4]}[$i]}=@infor2;
                                         } 
                            } 
           }
 close F2;  

}                            ##### end 3

############## print out
my @temp=keys %infor;
for(my $i=0;$i<@temp;$i++){
 for(my $j=0;$j<@{$infor{$temp[$i]}};$j++){
  print "@{$infor{$temp[$i]}[$j]}\n";
                                          }
}
