#!/usr/bin/perl -w
#use strict;


use Getopt::Std;
use vars qw/ $opt_f $opt_o /;
use List::Util qw(first);
my ($infile, $infile2, $barcode_length, $outprefix, $outfile2, $length, $name);

getopts('f:o:');  
&var_check();


open (IN, $infile) or die "Cannot open $infile: $!\n\n";
open (OUTP1, ">$outprefix\_P1.txt");



while (my $line = <IN>) {
  chomp $line;
    if($line =~ />/){
        my @header = split(/_/,$line);
        print OUTP1 "SNP_$header[-1]\tVVV\t11111111\t";
    } else{ 
            my $revline = reverse($line);
            my @SNPline = split(//,$line);
            my @revSNPline = split(//,$revline);
            my $linesize = length($line);
            my( $i ) = grep { evaluate($SNPline[$_]) } 0 .. $#SNPline;
            my( $h ) = grep { evaluate($revSNPline[$_]) } 0 .. $#revSNPline;
            my $fivePr_seq = substr($line,($i-20),15);
            my $threePr_seq = substr($line,(($linesize-$h)+5),15);
            my $fivep_end = ($i-20);
            my $threep_end = ($linesize - ($h-20));
            my $orig_size = ($threep_end-$fivep_end);
            my $orig_seq = substr($line,$fivep_end,$orig_size);
            
            print OUTP1 "$fivePr_seq\t$threePr_seq\t$orig_seq\n";
            #print "$fivep_end\t$threep_end\t$orig_size\n";
            #print "$orig_seq\n";
    }
}
close OUTP1;



sub evaluate {
    my($dna) = @_;
        if($dna =~ /[^ACGTacgt]/){
            return $dna;
        }
    
}

sub var_check {
        if ($opt_f) {
                $infile = $opt_f;
        } else {
                &var_error();
        }
            
        if ($opt_o) {
                $outprefix = $opt_o;
        } else {
                &var_error();
        } 
        
        
}

sub var_error {
        print "\n\n";
    
        print " REQUIRED:\n";
        print " -f     Input sequence.txt filename #1\n";
  
        print " -o     Output filename #1\n";
     
        print "\n\n";
        exit 0;
}


