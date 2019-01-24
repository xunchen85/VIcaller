#! usr/bin/perl -w
use strict;
die "perl $0 <lst><fa>\n" unless  @ARGV==2;
my $lst=$ARGV[0];
my $fa=$ARGV[1];

### read the gi
open IN,"$lst";
my %hash;
$hash{$lst}=1;
close IN;

### read the fasta
$fa=~/gz$/?(open IN,"gzip -cd $fa|"||die):(open IN,$fa||die);
$/=">";<IN>;$/="\n";
my %out;

while(<IN>){
 my $info=$1 if(/^(\S+)/);
 $/=">";
 my $seq=<IN>;
 $/="\n";
 $seq=~s/>|\r|\*//g;
 print ">$info\n$seq\n" if(exists $hash{$info} && ! exists $out{$info});
 $out{$info}=1;
}
close IN;

