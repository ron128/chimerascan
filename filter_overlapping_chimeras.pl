#!/lab/bin/perl
use Getopt::Std;
use Data::Dumper;
use strict;

my $usage = "Usage : filter_overlapping_chimeras.pl

The purpose of this file is to convert any eland multi-mapping alignment file into a BEDPE file
usable by BEDTools. The BEDPE format is as follows:
 
   1. chrom1 (chrom for the first end of the pair)
   2. start1 (start for the first end of the pair)
   3. end1 (end for the first end of the pair)
   4. chrom2 (chrom for the second end of the pair)
   5. start2 (start for the second end of the pair)
   6. end2 (end for the second end of the pair)
   7. name (custom name for the pair)
   8. score (custom score for the pair)
   9. strand1 (orientation for the first end of the pair. +/-)
  10. strand2 (orientation for the second end of the pair. +/-)
  11. any number of user-defined fields. Ignored, but carried-through by BEDTools 

The input for this script is:
-i Infile    [ mctp_30TUEAAXX_3_BEDPE_Gene.txt ]
-o Outfile   [ mctp_30TUEAAXX_3_filtered_genes_BEDPE.txt ]

\n\n";
########################################################################################
my %option=();
getopts("o:i:", \%option);
unless (defined $option{o}){die $usage}; my $outfile = $option{o};
unless (defined $option{i}){die $usage}; my $infile = $option{i};
#########################################################################################
#
# STEP 1: Import directories containing lanes to be analyzed
#
my $BEDPE_ISECT = $infile;

# Print out 'chimeric' mate pairs to a BEDPE file
my $OUT = $outfile;
open OUT, ">$OUT" or die "Can't open file $OUT";

my %FILTER;
open(BEDPE_ISECT, $BEDPE_ISECT);
while(<BEDPE_ISECT>){
    chomp;
    my(@v)=split(/\t/);
    my($c1) = $v[0];
    my($s1) = $v[1];
    my($e1) = $v[2];
    my($c2) = $v[3];
    my($s2) = $v[4];
    my($e2) = $v[5];
    my($mp) = $v[6];
    my($str1) = $v[8];
    my($str2) = $v[9];
    my($ch) = $v[12];
    my($sh) = $v[13];
    my($eh) = $v[14];

    
    if($c1 eq $ch && $c2 eq $ch){
	if( (($s1 >= $sh && $s1 <= $eh) || ($e1 >= $sh && $e1 <= $eh)) && (($s2 >= $sh && $s2 <= $eh) || ($e2 >= $sh && $e2 <= $eh))){
	    # Print to BAM file
	    $FILTER{$mp}[0]++;
	}
    }
}
close(BEDPE_ISECT);


#
# Now go through the file parsing only those that matter
#
my %CHIMERAS;
open(BEDPE_ISECT, $BEDPE_ISECT);
while(<BEDPE_ISECT>){
    #chomp;
    my(@v)=split(/\t/);
    my($c1) = $v[0];
    my($s1) = $v[1];
    my($e1) = $v[2];
    my($c2) = $v[3];
    my($s2) = $v[4];
    my($e2) = $v[5];
    my($mp) = $v[6];
    my($score) = $v[7];
    my($str1) = $v[8];
    my($str2) = $v[9];
    my($ch) = $v[12];
    my($sh) = $v[13];
    my($eh) = $v[14];
    my($ucsc) = $v[15];

    if($FILTER{$v[6]}[0] eq '' && $score eq '1'){
	print OUT;
	
	# Print BEDPE of candidates
#	foreach my $v (@v){
#	    print OUT "$v\t";
#	}
#	print OUT "\n";
    }
}
close(BEDPE_ISECT);

###################################################################################
#
# Christopher Maher
# chrmaher@med.mich.edu
# Last modified : 03/27/2008
#
###################################################################################
