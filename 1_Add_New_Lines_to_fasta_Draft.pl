#!/usr/bin/perl -w
use strict;


my $usage = "Usage: 1_Add_New_Lines_to_.fasta_Draft.pl <input file with returns> <fasta header> <out file name> <line length>\n\n";
my $fastq = shift or die $usage;
my $name = shift or die $usage;
my $out = shift or die $usage;
my $line_len =shift or die $usage;

my %genus;

open (DAT, $fastq) or die $!;
open (OUT, "> $out") or die $!;
my %count;
while (my $line = <DAT>) {
 chomp $line;
 if ($line =~ s/@/>/) {
    print OUT "$line\n";
    }  else {
   
    my $m = $line;
       $m =~ s/(.{1,$line_len})/$1\n/gs;
            print OUT "$m";


}
   
}

print OUT "\n";

close DAT;
close OUT;
