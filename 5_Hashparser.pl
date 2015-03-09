
#!/usr/bin/perl
use strict;
use warnings;


use Getopt::Std;
use vars qw/ $opt_f $opt_m $opt_o /;

my ($infile1, $infile2, $name, $min_hits, $v_conserved, $v_conserved_2, $r_conserved_2, $r_conserved, $outprefix1, $outprefix2, $a_header, $a_seq, $a_header2, $a_qual,
    $adapter_length, $seq_frag, $adapter_pos);



getopts('f:m:o:');
&var_check();

open (IN1, $infile1) or die "Cannot open $infile1: $!\n\n";





open (OUTP1, ">$outprefix1");


while(<IN1>){
    my $line = $_;
    chomp $line;
    my @SNPline = split('\t',$line);
    my @sample = split('\s',$SNPline[0]);
    
    
    
    if (($SNPline[-5] > $min_hits) && (exists($sample[1]))){
        print OUTP1 "$sample[0]\t$sample[1]\t$SNPline[-5]\t$SNPline[-1]\t$SNPline[-4]\n";
        
        
    }
        elsif((!exists($sample[1]) && ($SNPline[-5] > $min_hits))){
            print OUTP1 "$sample[0]\tNothing\t$SNPline[-1]\t$SNPline[-3]\t$SNPline[-2]\n";
    } else{
            
        }
    # print "$sample[0]\t$sample[1]\t$SNPline[-2]\n"; 
    #print OUTP1 "$SNPline[0]\t$SNPline[-2]\n";
   # if ($SNPline[-2] < $min_hits){
   # } else{
   #     print "$sample[0]\t$sample[1]\t$SNPline[-2]\n";
   # }
   # next;
}


            
           
print "\n";


close OUTP1;









sub var_check {
        if ($opt_f) {
                $infile1 = $opt_f;
        } else {
                &var_error();
        }
        
       	if ($opt_m) {
                $min_hits = $opt_m;
        } else {
                &var_error();
        }
		if ($opt_o) {
                $outprefix1 = $opt_o;
        } else {
                &var_error();
        }
        
		
        
}

sub var_error {
        print "\n\n";
        
        print " You did not provide enough information.\n\n";
     
        print " REQUIRED:\n";
        print " -f     	Hash infile\n";
		print " -m		minimum number of hits to keep\n";
        print " -o     	Output filename prefix i.e. MavA5\n";
        print "\n\n";
        exit 0;
}
