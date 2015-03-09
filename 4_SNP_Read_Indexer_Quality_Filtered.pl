#!/usr/bin/perl
use strict;
use warnings;                                                                                       
#####################################################################################################
#							#####The SNP Read Grabber#####								            #
#                                                                                                   #
#   input1 = illumina.fastq                                                                         #
#   input2 = SNP file (tab delimited ## SNP120	5554678	R	GAGAAGGATTGGCCAGGA	GTCGAATCTTCGGTTCGC  #
#	output1-n = File for every SNP site containing a hash of the variation	                        #
#                                                                                                   #
#####################################################################################################

use Getopt::Std;
use vars qw/ $opt_f $opt_F $opt_o $opt_q /;


my ($infile_Reads, $infile_SNPs, $outprefix, $seq_frag, $adapter_pos, $qual_thresh);

getopts('f:F:o:q:');
&var_check();


open (IN2, $infile_SNPs) or die "Cannot open $infile_SNPs: $!\n\n";
open (OUT, ">$outprefix");


while(my $line = <IN2>){
    chomp $line;
    
	my @SNPline = split('\t', $line);
    
	my $fivePr_seq = $SNPline[3];                    # Creats variable from  5' flanking sequence from SNP file
        my $threePr_seq = $SNPline[4];                   # Creats variable from  3' flanking sequence from SNP file
	                                            
	my $revcomp_fivePr = revcom($fivePr_seq);		 # Reverse complements flanking sequences into variable
	my $revcomp_threePr = revcom($threePr_seq);      # 					"
	       
	  
	                                                
    my $adapter_length = length($fivePr_seq);        # Puts length of Flanking seqs into variables
         
    my $variable_length;
	my %SNP;
	open (IN1, $infile_Reads) or die "Cannot open $infile_Reads: $!\n\n";
	
    while(<IN1>) {
        my $header = $_;
        my $seq = <IN1>;
        my $header2 = <IN1>;       
        my $qual = <IN1>;
	#print "header2=$header2,seq=$seq,header=$header,qual=$qual\n";
        chomp $seq;
        chomp $header2;
        chomp $header;
        chomp $qual;
	my $seq_length = length($seq);
	
	
	
        if (($seq =~ /$fivePr_seq/i) && ($seq =~ /$threePr_seq/i)) {        # Pulls out reads containing the two flanking seqs avoding grabbing same seq twice
            my $conserved_position_1 = substr($seq, (index ($seq, $fivePr_seq) + $adapter_length));   # variable of all characters downstream of 5' seq
            my $conserved_position_2 = substr($conserved_position_1, 0, index($conserved_position_1, $threePr_seq));
			$variable_length = length($conserved_position_2);
			
			my $qual_frag1 = length($conserved_position_1);
			
			my @phred = split(//,(substr($qual, ($seq_length - $qual_frag1),$variable_length)));
			@phred = map(ord, @phred);
			my $score = 0;
			foreach(@phred){
				if($_ >= $qual_thresh){
					$score++;
				}
			}
			if($score == $variable_length){
					#print "$conserved_position_2\n";
                if (exists($SNP{$conserved_position_2})) {
                    $SNP{$conserved_position_2}++;
                }  else{
                    $SNP{$conserved_position_2} = 1;
      } 
                    
			}
				else{				
				}
				
			
	} elsif (($seq =~ /$threePr_seq/i) && ($seq !~ /$fivePr_seq/i)) { 
            my $conserved_position_1 = substr($seq, 0, index($seq, $threePr_seq)); # variable of all characters upstream of 3' seq
			my $conserved_position_2 = substr($conserved_position_1, (index($seq, $threePr_seq) - ($variable_length)));
			
			my $qual_frag1 = length($conserved_position_1);
			my @phred = split(//,(substr($qual, ($seq_length - $qual_frag1),$variable_length)));
			@phred = map(ord, @phred);
			my $score = 0;
			foreach(@phred){
				if($_ >= $qual_thresh){
					$score++;
				}
			}
				if ((length($conserved_position_2) >= $variable_length) && ($score == $variable_length)){
					#print "$conserved_position_2\n";
				if (exists($SNP{$conserved_position_2})) {
                    $SNP{$conserved_position_2}++;
                }  else{
                    $SNP{$conserved_position_2} = 1;
                    } 
				}
					else{				
					}			
        } elsif (($seq !~ /$threePr_seq/i) && ($seq =~ /$fivePr_seq/i)) {
            my $conserved_position_1 = substr($seq, (index ($seq, $fivePr_seq) + $adapter_length)); # variable of all characters downstream of 5' seq
            my $conserved_position_2 = substr($conserved_position_1, 0,($variable_length));
			
			my $qual_frag1 = length($conserved_position_1);
			my @phred = split(//,(substr($qual, ($seq_length - $qual_frag1),$variable_length)));
			@phred = map(ord, @phred);
			my $score = 0;
			foreach(@phred){
				if($_ >= $qual_thresh){
					$score++;
				}
			}
				if ((length($conserved_position_2) >= $variable_length) && ($score == $variable_length)){
					#print "$conserved_position_2\n";
				if (exists($SNP{$conserved_position_2})) {
                    $SNP{$conserved_position_2}++;
                }  else{
                    $SNP{$conserved_position_2} = 1;
                    } 				
				}
					else{				
					}			
        } elsif (($seq =~ /$revcomp_fivePr/i) && ($seq !~ /$revcomp_threePr/i)) {
            my $rev_seq = revcom($seq);
            my $conserved_position_1 = substr($rev_seq, (index ($rev_seq, $fivePr_seq) + $adapter_length));
			my $conserved_position_2 = substr($conserved_position_1, 0,($variable_length));
			
			my $qual_frag1 = length($conserved_position_1);
			my @phred = split(//,(substr($qual, ($seq_length - $qual_frag1),$variable_length)));
			@phred = map(ord, @phred);
			my $score = 0;
			foreach(@phred){
				if($_ >= $qual_thresh){
					$score++;
				}
			}
				if ((length($conserved_position_2) >= $variable_length) && ($score == $variable_length)){
					#print "$conserved_position_2\n";
				if (exists($SNP{$conserved_position_2})) {
                    $SNP{$conserved_position_2}++;
                }  else{
                    $SNP{$conserved_position_2} = 1;
      } 				
				}
					else{				
					}					
        } elsif (($seq =~ /$revcomp_threePr/i) && ($seq !~ /$revcomp_fivePr/i)) { 
            my $rev_seq = revcom($seq);
			my $conserved_position_1 = substr($rev_seq, 0, index($rev_seq, $threePr_seq));
			my $conserved_position_2 = substr($conserved_position_1, (index($rev_seq, $threePr_seq) - ($variable_length)));
			
			my $qual_frag1 = length($conserved_position_1);
			my @phred = split(//,(substr($qual, ($seq_length - $qual_frag1),$variable_length)));
			@phred = map(ord, @phred);
			my $score = 0;
			foreach(@phred){
				if($_ >= $qual_thresh){
					$score++;
				}
			}
				if ((length($conserved_position_2) >= $variable_length) && ($score == $variable_length)){
					#print "$conserved_position_2\n";
                    if (exists($SNP{$conserved_position_2})) {
                        $SNP{$conserved_position_2}++;
                    }  else{
                        $SNP{$conserved_position_2} = 1;
                        } 
                    }
					else{						
					}
			
        } elsif (($seq =~ /$revcomp_threePr/i) && ($seq =~ /$revcomp_fivePr/i)) {
            my $rev_seq = revcom($seq);
            my $conserved_position_1 = substr($rev_seq, (index ($rev_seq, $fivePr_seq) + $adapter_length));
			my $conserved_position_2 = substr($conserved_position_1, 0, index($conserved_position_1, $threePr_seq));
			
			my $qual_frag1 = length($conserved_position_1);
			my @phred = split(//,(substr($qual, ($seq_length - $qual_frag1),$variable_length)));
			@phred = map(ord, @phred);
			my $score = 0;
			foreach(@phred){
				if($_ >= $qual_thresh){
					$score++;
				}
			}
				if ((length($conserved_position_2) >= $variable_length) && ($score == $variable_length)){
					#print "$conserved_position_2\n";
                    if (exists($SNP{$conserved_position_2})) {
                        $SNP{$conserved_position_2}++;
                    }  else{
                        $SNP{$conserved_position_2} = 1;
                        } 
				}
					else{
					}
                    
        }
            
        else {
                    
			next;
        }
		
	}
					for my $key ( keys %SNP ) {
						my $value = $SNP{$key};
						print OUT  "$SNPline[0] $key\t$value\t$SNPline[3]$key$SNPline[4]\t$SNPline[3]\t$SNPline[4]\t$SNPline[5]\n";
					}
                   
             
} 
close IN2;
close OUT;

sub revcom {
    my($dna) = @_;
    my $revcom = reverse $dna;
    $revcom =~ tr/ACGTacgt/TGCAtgca/;
    return $revcom;
}





sub var_check {
        if ($opt_f) {
                $infile_Reads = $opt_f;
        } else {
                &var_error();
        }
          if ($opt_F) {
                $infile_SNPs = $opt_F;
        } else {
                &var_error();
        }
		    if ($opt_o) {
                $outprefix = $opt_o;
        } else {
                &var_error();
        }
          if ($opt_q) {
                $qual_thresh = $opt_q;
        } else {
                &var_error();
        }
		  
}

sub var_error {
        print "\n\n";
        
        print " You did not provide enough information.\n\n";
     
        print " REQUIRED:\n";
        print " -f     Read file neame\n";
		print " -F     SNP file name\n";
        print " -o     Output filename prefix\n";
		print " -q     Quality Score Threshold\n";
        print "\n\n";
        exit 0;
}
