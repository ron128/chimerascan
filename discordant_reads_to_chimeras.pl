#!/usr/bin/perl

use Getopt::Std;
use Bio::DB::Sam;
use strict;

my $usage = "Usage : discordant_reads_to_chimeras.pl 

-b | Bam input file
-o | Output directory

Input is a discordant set of mate pairs; therefore neither are labeled as 'is pair mapped'

Sample Input
/exds/users/mkiyer/projects/fusiondiscovery/example3/discordant_reads.bam (2 x 53)

Potential updates:
=====================
- Might want to filter mitochondrial before choosing best hit


References
# http://search.cpan.org/~lds/Bio-SamTools/lib/Bio/DB/Bam/Alignment.pm
# http://kobesearch.cpan.org/htdocs/Bio-SamTools/Bio/DB/Sam.pm.html#features_sam_gt_features_options

\n\n";

################################ VARIABLES AND SETTINGS #############################
my %option=();
getopts("b:f:o:", \%option);
unless (defined $option{b}){die $usage}; my $bamfile = $option{b};
#unless (defined $option{f}){die $usage}; my $flowid = $option{f};
unless (defined $option{o}){die $usage}; my $outdir = $option{o};
#####################################################################################
#
# Configure bam file
#
my $bam = Bio::DB::Sam->new(-bam =>$bamfile);


#
# Output files
#

# Print out 'Single Best Split' mate pairs to a BEDPE file
my $OUTS = $outdir.'/split_chimeras.bedpe.txt';
open OUTS, ">$OUTS" or die "Can't open file $OUTS";

# Print out 'Multi-mapping split' mate pairs to a BEDPE file
#my $OUTMM = $outdir.'/multimap_split_chimeras.bedpe.txt';
#open OUTMM, ">$OUTMM" or die "Can't open file $OUTMM";

# Print out 'One-mapper' mate pairs to a BEDPE file
my $OUTONE = $outdir.'/onemappers_chimeras.bedpe.txt';
open OUTONE, ">$OUTONE" or die "Can't open file $OUTONE";

# Print out 'Non-mapper' mate pairs to a BEDPE file
my $OUTNON = $outdir.'/nonmappers_chimeras.bedpe.txt';
open OUTNON, ">$OUTNON" or die "Can't open file $OUTNON";




#
# Establish counters
#
my $total_align_attempt = '0';
my $both_non_mapping = '0';
my $single_non_mapping = '0';
my $split_mates = '0';
my $current_mate_id = '';
my (%READ1,%READ2);



#
# Process mate alignments
#
my $iterator = $bam->features(-iterator=>1); 
while (my $a = $iterator->next_seq) {
    my $read_name = $a->query->name;
   
    # Confirm if processing new read
    if($current_mate_id ne '' && $read_name ne $current_mate_id){ 
	process_alignments(\%READ1,\%READ2,$read_name,$current_mate_id); 
	%READ1 = ();
	%READ2 = ();
    }
    $current_mate_id = $read_name;

    my $seqid  = $a->seq_id;
    my $start  = $a->start;
    my $end    = $a->end;
    my $strand = $a->strand;
#    my $ref_dna= $a->dna; => Gives error
    my $query_start  = $a->query->start;
    my $query_end    = $a->query->end;
    my $query_strand = $a->query->strand;
    my $query_dna    = $a->query->dna;
    my $bitflag   = $a->flag;
    my $cigar     = $a->cigar_str;
    my @scores    = $a->qscore;     # per-base quality scores
    my $match_qual= $a->qual;       # quality of the match
    my $is_paired = $a->paired;     # Return true if the aligned read is part of a mate/read pair (regardless of whether the mate mapped).
    my $is_proper = $a->proper_pair;   # Return true if the aligned read is part of a mate/read pair and both partners mapped to the reference sequence.
    my $is_unmapped = $a->unmapped;    # Return true if the read failed to align.
    my $mate_is_unmapped = $a->munmapped; # Return true if the read's mate failed to align.
    my $reversed = $a->reversed;       # Return true if the aligned read was reverse complemented prior to aligning.
    my $mate_reversed = $a->mreversed; # Return true if the aligned read's mate was reverse complemented prior to aligning.
#    my $mseqid = $a->mate_seq_id;      # Return the seqid of the mate.
#    my $mstart = $a->mate_start;       # For paired reads, return the start of the mate's alignment in reference sequence coordinates.
#    my $mend = $a->mate_end; # For paired reads, return the end position of the mate's alignment. in reference sequence coordinates.
                             # -item $len = $align->mate_len
                             # For mate-pairs, retrieve the length of the mate's alignment on the reference sequence.
    my $insert = $a->isize;
    my $str = $a->aux;
    my @tags = $a->get_all_tags;
    my @values = $a->get_tag_values('FLAGS');
    my $is_duplicate_true = $a->has_tag('FIRST_MATE');
    
    # Determine mate
    my $first_mate = '0';
    my $second_mate = '0';
    foreach my $v (@values){ if($v =~ /FIRST/){ $first_mate++; }else{ $second_mate++; } }

    # Extract alignment info (mismatches, alignments)
    my($mismatches,$mappings) = extract_alignment_info($str);

    # Paired (Always paired in this format)
    # Is_unmapped (1 - not mapped; 0 - is mapped)
    if($is_paired eq '1' && $is_unmapped eq '0'){
	
	# Filter PhiX / HSU
	if($seqid =~ /phiX/ || $seqid =~ /HSU/){
	    next;
	}

	# Push alignment
	my $id = $seqid.'_'.$start.'_'.$end.'_'.$strand.'_'.$insert.'_'.$mismatches.'_'.$mappings.'_'.$query_dna;
	if($first_mate eq '1'){
	    push @{$READ1{$read_name}}, $id;
	}else{
	    push @{$READ2{$read_name}}, $id;
	}

	

    }elsif($is_paired eq '1' && $is_unmapped eq '1'){
	
	# Push alignment
	my $id = 'None_None_None_None_'.$insert.'_None_'.$mappings;
	
	if($first_mate eq '1'){
	    push @{$READ1{$read_name}}, $id;
	}else{
	    push @{$READ2{$read_name}}, $id;
	}
    }else{
	print "ERROR: Is not paired\n";
    }
}

close(OUTSB);
close(OUTMM);
close(OUTONE);
close(OUTNON);

###################################################################################
sub process_alignments{
    my($R1_href,$R2_href,$read_name,$current_mate_id) = @_;
 
    # PRINT TO BEDPE FILE FORMAT
    # 1. chrom1 (chrom for the first end of the pair)
    # 2. start1 (start for the first end of the pair)
    # 3. end1 (end for the first end of the pair)
    # 4. chrom2 (chrom for the second end of the pair)
    # 5. start2 (start for the second end of the pair)
    # 6. end2 (end for the second end of the pair)
    # 7. name (custom name for the pair)
    # 8. score (custom score for the pair)
    # 9. strand1 (orientation for the first end of the pair. +/-)
    # 10. strand2 (orientation for the second end of the pair. +/-)
    # 11. any number of user-defined fields. Ignored, but carried-through by BEDTools 

    if($$R1_href{$current_mate_id}[0] ne '' && $$R2_href{$current_mate_id}[0] ne ''){
	
	# NON-Mappings
	if($$R1_href{$current_mate_id}[0] =~ /None/ && $$R2_href{$current_mate_id}[0] =~ /None/){
	    my($seqid1,$start1,$end1,$strand1,$insert1,$mismatches1,$mappings1,$seq1) = split(/\_/,$$R1_href{$current_mate_id}[0]);
	    my($seqid2,$start2,$end2,$strand2,$insert2,$mismatches2,$mappings2,$seq2) = split(/\_/,$$R2_href{$current_mate_id}[0]);
	    print OUTNON "$seqid1\t$start1\t$end1\t$seqid2\t$start2\t$end2\t$current_mate_id\t1\t$strand1\t$strand2\tSPLIT;$seq1;$seq2\n";
	    return;
	}


	# ONE-Mappers
	if($$R1_href{$current_mate_id}[0] =~ /None/ || $$R2_href{$current_mate_id}[0] =~ /None/){
	    my($seqid1,$start1,$end1,$strand1,$insert1,$mismatches1,$mappings1,$seq1) = split(/\_/,$$R1_href{$current_mate_id}[0]);
	    my($seqid2,$start2,$end2,$strand2,$insert2,$mismatches2,$mappings2,$seq2) = split(/\_/,$$R2_href{$current_mate_id}[0]);
	    print OUTONE "$seqid1\t$start1\t$end1\t$seqid2\t$start2\t$end2\t$current_mate_id\t1\t$strand1\t$strand2\tSPLIT;$seq1;$seq2\n";
	    return;
	}


	# Split mate pairs 	

	my($best1,$best1_count) = select_best_mate($R1_href,$current_mate_id);
	my($best2,$best2_count) = select_best_mate($R2_href,$current_mate_id);

	if($best1_count eq '1' && $best2_count eq '1'){
	    
	    # Iterate through mappings to chose best pairings
	    for(my $x = '0'; $x < @{$$R1_href{$current_mate_id}}; $x++){
		for(my $y = '0'; $y < @{$$R2_href{$current_mate_id}}; $y++){
		
		    if($$R1_href{$current_mate_id}[$x] =~ /chrM/ || $$R2_href{$current_mate_id}[$y] =~ /chrM/){ 
			last;
		    }
		    
		    if($$R1_href{$current_mate_id}[$x] =~ /illumina/ || $$R2_href{$current_mate_id}[$y] =~ /illumina/){ 
			last;
		    }
		    
		    my($seqid1,$start1,$end1,$strand1,$insert1,$mismatches1,$mappings1,$seq1) = split(/\_/,$$R1_href{$current_mate_id}[$x]);
		    my($seqid2,$start2,$end2,$strand2,$insert2,$mismatches2,$mappings2,$seq2) = split(/\_/,$$R2_href{$current_mate_id}[$y]);
		    
		    #
		    # Since best equals 1, need to print only those that have the best number of mismatches
		    #
		    if($mismatches1 eq $best1 && $mismatches2 eq $best2){
			print OUTS "$seqid1\t$start1\t$end1\t$seqid2\t$start2\t$end2\t$current_mate_id\t1\t$strand1\t$strand2\tSINGLE;$seq1;$seq2\n";
		    }
		}
	    }
	}else{

	    #
	    # Figure out what % are multi-mappers
	    #
	    
	    # Iterate through mappings to chose best pairings
	    for(my $x = '0'; $x < @{$$R1_href{$current_mate_id}}; $x++){
		for(my $y = '0'; $y < @{$$R2_href{$current_mate_id}}; $y++){
		    
		    if($$R1_href{$current_mate_id}[$x] =~ /chrM/ || $$R2_href{$current_mate_id}[$y] =~ /chrM/){ 
			last;
		    }
		    
		    if($$R1_href{$current_mate_id}[$x] =~ /illumina/ || $$R2_href{$current_mate_id}[$y] =~ /illumina/){ 
			last;
		    }
		    
		    my($seqid1,$start1,$end1,$strand1,$insert1,$mismatches1,$mappings1,$seq1) = split(/\_/,$$R1_href{$current_mate_id}[$x]);
		    my($seqid2,$start2,$end2,$strand2,$insert2,$mismatches2,$mappings2,$seq2) = split(/\_/,$$R2_href{$current_mate_id}[$y]);
		    
		    #
		    # Print pairs having the lowest (best) number of mismatches
		    #
		    if($mismatches1 eq $best1 && $mismatches2 eq $best2){
			print OUTS "$seqid1\t$start1\t$end1\t$seqid2\t$start2\t$end2\t$current_mate_id\t1\t$strand1\t$strand2\tMULTI;$seq1;$seq2\n";
		    }
		}
	    }
	}
    }
}

sub select_best_mate{
    my($R_href,$current_mate_id) = @_;

    my $best = ''; # Lower better
    my $best_count = '0'; # Determine if its a single best

    # Iterate through each mapping
    for(my $x = '0'; $x < @{$$R_href{$current_mate_id}}; $x++){
	if($$R_href{$current_mate_id}[$x] =~ /chrM/){ 
	    return('fail');
	}

	my($seqid,$start,$end,$strand,$insert,$mismatches,$mappings) = split(/\_/,$$R_href{$current_mate_id}[$x]);
	if($best eq '' || $mismatches < $best){
	    $best = $mismatches;
	    $best_count = '1';
	}elsif($mismatches eq $best){
	    $best_count++;
	}
#	print "$current_mate_id\t$seqid\t$mismatches\t$best\t$best_count\n";
    }    
    return($best,$best_count);
}

sub extract_alignment_info{
    my($str) = @_;
    
    my($mismatches,$mappings);
    my(@mappings) = split(/\s+/,$str);
    foreach my $map (@mappings){ 
	if($map =~ /^NH/){
	    my($NH,$i,$maps) = split(/\:/,$map);
	    $mappings = $maps;
	}
	if($map =~ /^NM/){
	    my($NM,$i,$mis) = split(/\:/,$map);
	    $mismatches = $mis;
	}
    }
    
    return($mismatches,$mappings);
}

###################################################################################
#
# Christopher Maher
# chrmaher@med.mich.edu
# Last modified : 09/21/2007
#
###################################################################################
