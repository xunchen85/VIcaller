#!usr/bin/perl
use strict;
open FASTA2,"$ARGV[0]";
my $line="";
my @line=();
my $read_name="";
my $read_seq="";
my %fasta2=();

################## read fasta2 file;
while(<FASTA2>){
 $line=$_;
 chmod($line);
 $line=~s/\s+$//;
 if($line=~ "^>"){
  $read_name=$line;
  $fasta2{$read_name}="";
                 }
 else{
  $read_seq=$line;
  if($fasta2{$read_name}){$fasta2{$read_name}.=$read_seq;}
  else{$fasta2{$read_name}=$read_seq;}
     }
}

my @names=keys %fasta2;
@names=sort @names;
for(my $i=0;$i<@names;$i++){
 my $count = () = $fasta2{$names[$i]} =~ /C|A|T|G|a|t|c|g/g;
 if($count>=20){print "$names[$i]\n$fasta2{$names[$i]}\n";}
}
