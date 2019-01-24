#!usr/bin/perl
use strict;
my $group=$ARGV[1];

############# Extract read seq
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
#   print "$read_name\t$temp04\n";
}
close F1;

open F1,"$ARGV[0]_alinged_sf.fuq"||di("di");
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

################# repeat;
my %repeats1=();
my %repeats2=();
open F3,"$ARGV[0]_alinged_both.repeat2";
while(<F3>){
 my @line=split;
 my $read_name=$line[0];
 $read_name=~s/\/[12]\s+$//;
 $read_name=~s/^@//;
 $read_name=~s/\/1//;
 $read_name=~s/\/2//;
 if($line[0] =~ "/1"){@{$repeats1{$read_name}}=@line;}
 elsif($line[0] =~ "/2"){@{$repeats2{$read_name}}=@line;}
}

################# read the file containing full information of each chimeric and split read;
my @infor1=();
my @infor2=();
my %infor=();
open F1,"$ARGV[0].information";
while(<F1>){
 @infor1=split;
 my $name1=join(" ",@infor1);
 @{$infor{$name1}}=@infor1;
 if(exists($file{$infor1[7]})){push(@{$infor{$name1}},$file{$infor1[7]});} else {push(@{$infor{$name1}},"NA");}
 if(exists($file2{$infor1[7]})){push(@{$infor{$name1}},$file2{$infor1[7]});} else {push(@{$infor{$name1}},"NA");}
 if(exists($repeats1{$infor1[7]})){push(@{$infor{$name1}},@{$repeats1{$infor1[7]}}[1..4]);} else {push(@{$infor{$name1}},("NA","NA","NA","NA"));}
 if(exists($repeats2{$infor1[7]})){push(@{$infor{$name1}},@{$repeats2{$infor1[7]}}[1..4]);} else {push(@{$infor{$name1}},("NA","NA","NA","NA"));}
}
       
######################## target default;
my $line="";
my @line=();
my %f1=();
open F1,"$ARGV[0]_aligned_target$group";
while(<F1>){
 @line=split;
 my $name2=join(" ",@line[0..13]);
 @{$f1{$name2}}=@line;
}
my @temp1=keys %infor;
for(my $i=0;$i<@temp1;$i++){
 if(exists($f1{$temp1[$i]})){push(@{$infor{$temp1[$i]}},@{$f1{$temp1[$i]}}[14..58]);}
 else{ my @empty= (0) x 45; push(@{$infor{$temp1[$i]}},@empty);}
}
close F1;

###################### vicaller;
my %f1=();
open F1,"$ARGV[0]_aligned_vicaller$group";
while(<F1>){
 @line=split;
 my $name2=join(" ",@line[0..13]);
 @{$f1{$name2}}=@line;
} 
my @temp1=keys %infor;
for(my $i=0;$i<@temp1;$i++){
 if(exists($f1{$temp1[$i]})){push(@{$infor{$temp1[$i]}},@{$f1{$temp1[$i]}}[14..58]);}
 else{ my @empty= (0) x 45; push(@{$infor{$temp1[$i]}},@empty);}
}
close F1;

###################### human;
my %f1=();
open F1,"$ARGV[0]_aligned_human$group";
while(<F1>){
 @line=split;
 my $name2=join(" ",@line[0..13]);
 @{$f1{$name2}}=@line;
}
my @temp1=keys %infor;
for(my $i=0;$i<@temp1;$i++){
 if(exists($f1{$temp1[$i]})){push(@{$infor{$temp1[$i]}},@{$f1{$temp1[$i]}}[14..58]);}
 else{ my @empty= (0) x 45; push(@{$infor{$temp1[$i]}},@empty);}
}
close F1;

#####################
my @temp1=keys %infor;
for(my $i=0;$i<@temp1;$i++){
 print "@{$infor{$temp1[$i]}}\n";
}
