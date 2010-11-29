#!/lab/bin/perl
use Getopt::Std;
use Data::Dumper;
use strict;

my $usage = "Usage : extract_chimera_candidates.pl

The input for this script is:
-o Outdir    [ /home/chrmaher/CHIMERASCAN_1.0/CANDIDATES/ ]
-i Input     [ /home/chrmaher/CHIMERASCAN_1.0/FILTERED4GENES_BEDPE/30TUEAAXX_s_3_filtered_genes_BEDPE.txt ]
-u UCSC      [ UCSC_genes_file.txt ]
# Passing filter
# Non-redundant mappings

\n\n";
########################################################################################
my %option=();
getopts("f:o:i:u:", \%option);
unless (defined $option{o}){die $usage}; my $outdir = $option{o};
unless (defined $option{i}){die $usage}; my $IN = $option{i};
unless (defined $option{u}){die $usage}; my $ucscfile = $option{u};
#########################################################################################
#
# STEP 1: Extract UCSC gene positions
#
my($UCSC_INFO,$UCSC2GENE,$UCSCID,$ORDER,$AORDER,$EXONS) = get_ucsc($ucscfile);

# Print out 'chimeric' mate pairs to a tab delimited file
my $OUT = $outdir.'filtered_chimeras.txt';
#print "OUT $OUT\n";
open OUT, ">$OUT" or die "Can't open file $OUT";

# Print out 'chimeric' mate pairs to a BEDPE file
my $OUTBEDPE = $outdir.'filtered_chimeras.bedpe.txt';
#print "OUTBEDPE $OUTBEDPE\n";
open OUTBEDPE, ">$OUTBEDPE" or die "Can't open file $OUTBEDPE";

# Print out 'chimeric' nominations into FASTA file or reads
my $READS = $outdir.'filtered_chimeras_reads.fa';
#print "READS $READS\n";
open READS, ">$READS" or die "Can't open file $READS";


#
# STEP 2: PARSE FILE
#
my(%SEQS5P,%SEQS3P,%POS,%READS5P_TMP,%READS3P_TMP,%READS5P,%READS3P);
#print "$IN\n";
open(IN, $IN);
while(<IN>){
    chomp;
    my(@v)=split(/\t/);
    my($c1) = $v[0];
    my($s1) = $v[1];
    my($e1) = $v[2];
    my($c2) = $v[3];
    my($s2) = $v[4];
    my($e2) = $v[5];
    my($mp) = $v[6];
    my($str1) = $v[8]; if($str1 eq 'F' || $str1 eq '1'){ $str1 = '+'; } if($str1 eq 'R' || $str1 eq '-1'){$str1 = '-'; }
    my($str2) = $v[9]; if($str2 eq 'F' || $str2 eq '1'){ $str2 = '+'; } if($str2 eq 'R' || $str2 eq '-1'){$str2 = '-'; }
    my($ch) = $v[11];
    my($sh) = $v[12];
    my($eh) = $v[13];
    my($ucsc) = $v[14];
    my($ucscstr) = $v[16];
    my($exons) = $v[20];
    my(@exon_lengths) = split(/\,/,$v[21]); 
    my(@exon_starts) = split(/\,/,$v[22]);

    # If first mate matches hit
    if($c1 eq $ch){
	
	# If first mate falls within boundaries of hit
	if( ($s1 >= $sh && $s1 <= $eh) || ($s1 >= $sh && $s1 <= $eh) ){
	    my $gene = $$UCSC2GENE{$ucsc}[0];
		print "$gene\n";

	    # Can determine if 5' or 3'
 	    my $partner;
	    if($str1 eq $ucscstr){
		$partner = '5p';
	    }else{
		$partner = '3p';
	    }
	    my($mate_position,$exon) = determine_position($partner,$exons,$ch,$sh,$eh,\@exon_lengths,\@exon_starts,$ucscstr,$s1,$e1,$ucsc);
	     if($mate_position eq ''){
		 next;
	     }else{
		 my $tmp = $gene.';'.$ucsc.';'.$partner.';'.$mate_position.';'.$ucscstr.';'.$exon;
		 my($match,$read1,$read2) = split(/\;/,$v[10]);
		 if($match eq 'SINGLE'){
		     #push @{$SEQS{$mp}},$read1,$read2;
		      push @{$SEQS5P{$mp}},$read1;
		      push @{$SEQS3P{$mp}},$read2;
		     if($partner eq '5p'){ push @{$READS5P{$mp}}, $tmp; }elsif($partner eq '3p'){  push @{$READS3P{$mp}}, $tmp; }
		 }
	     }
        # If second mate falls within boundaries of hit
	}elsif($c2 eq $ch){
	    if( ($s2 >= $sh && $s2 <= $eh) || ($s2 >= $sh && $s2 <= $eh) ){
		my $gene = $$UCSC2GENE{$ucsc}[0];
		# Can determine if 5' or 3'
		my $partner;
		if($str2 eq $ucscstr){
		    $partner = '5p';
		}else{
		    $partner = '3p';
		}
		my($mate_position,$exon) = determine_position($partner,$exons,$ch,$sh,$eh,\@exon_lengths,\@exon_starts,$ucscstr,$s2,$e2,$ucsc);
		 if($mate_position eq ''){
		     next;
		 }else{
		     my $tmp = $gene.';'.$ucsc.';'.$partner.';'.$mate_position.';'.$ucscstr.';'.$exon;
		     my($match,$read1,$read2) = split(/\;/,$v[10]);
		      if($match eq 'SINGLE'){
			  push @{$SEQS5P{$mp}},$read1;
			  push @{$SEQS3P{$mp}},$read2;
			  if($partner eq '5p'){ push @{$READS5P{$mp}}, $tmp; }elsif($partner eq '3p'){  push @{$READS3P{$mp}}, $tmp; }
		      }
		 }
	    }
	}    
    # If second mate matches hit	
    }elsif($c2 eq $ch){
	if( ($s2 >= $sh && $s2 <= $eh) || ($s2 >= $sh && $s2 <= $eh) ){
	    my $gene = $$UCSC2GENE{$ucsc}[0];
	    
	    # Can determine if 5' or 3'
	    my $partner;
	    if($str2 eq $ucscstr){
		$partner = '5p';
	    }else{
		$partner = '3p';
	    }


	    # Extract exon start and end
	    my($mate_position,$exon,$exon_start,$exon_end) = determine_position($partner,$exons,$ch,$sh,$eh,\@exon_lengths,\@exon_starts,$ucscstr,$s2,$e2,$ucsc);
	    if($mate_position eq ''){
		next;
	    }else{
		my $tmp = $gene.';'.$ucsc.';'.$partner.';'.$mate_position.';'.$ucscstr.';'.$exon.';'.$exon_start.';'.$exon_end;
		my($match,$read1,$read2) = split(/\;/,$v[10]);
		if($match eq 'SINGLE'){
		    #push @{$SEQS{$mp}},$read1,$read2;
		     push @{$SEQS5P{$mp}},$read1;
		     push @{$SEQS3P{$mp}},$read2;
		    if($partner eq '5p'){ push @{$READS5P{$mp}}, $tmp; }elsif($partner eq '3p'){  push @{$READS3P{$mp}}, $tmp; }
		}
	    }
	}
    }    
}
close(IN);

#my $FPEE = '2';
#my $FPucsc = 'uc010vof.1';
#my($relstart) = extract_rel_exon('5p',$FPEE,$FPucsc,$EXONS);
#print "$relstart\n";
#die;


#
#
# Filter each read
#
#
my(%FINAL,%FPpositions,%TPpositions,%TPexons,%FPexons,%MP,%FPexonstarts,%TPexonstarts,%FPexonends,%TPexonends);
my @mps = keys %SEQS5P;
foreach my $mp (@mps){

    #if($mp eq 'PATHBIO-SOLEXA2_30TUEAAXX:3:3:1614:1230'){
	#print "WOOOO\t$mp\n";
	#    }

    if($READS5P{$mp}[0] eq '' || $READS3P{$mp}[0] eq ''){
	next;
    }
    
    # Process 5p
    for(my $i = '0'; $i < @{$READS5P{$mp}}; $i++){

	for(my $j = '0'; $j < @{$READS3P{$mp}}; $j++){
	    
	    my(@x) = split(/\;/,$READS5P{$mp}[$i]); # $gene.';'.$ucsc.';'.$partner.';'.$mate_position.';'.$ucscstr.';'.$exon.';'.$exon_start.';'.$exon_end
	    my(@y) = split(/\;/,$READS3P{$mp}[$j]); # $gene.';'.$ucsc.';'.$partner.';'.$mate_position.';'.$ucscstr.';'.$exon.';'.$exon_start.';'.$exon_end
	    
	    #$UCSC{$x[1]}[0]++;
	    #$UCSC{$x[1]}[1] = $x[0]; # HUGO Gene ID 
	    #$UCSC{$x[1]}[2] = $x[3]; # Position
	    #$UCSC{$x[1]}[3] = $x[2]; # Partner
	    #$UCSC{$x[1]}[4] = $x[5]; # Exon
	
	    # Are these overlapping genes with different gene names?
	    my($fchr,$fstr,$fst,$fend) = split(/\;/,$$UCSCID{$x[1]}[0]);
	    my($tchr,$tstr,$tst,$tend) = split(/\;/,$$UCSCID{$y[1]}[0]);
		
	    if($tchr eq $fchr){
		if( ($tst >= $fst && $tst <= $fend) && ($tend >= $fst && $tend <= $fend) ){
		    next;
		}
		if( ($fst >= $tst && $fst <= $tend) && ($fend >= $tst && $fend <= $tend) ){
		    next;
		}		
	    }    
	    my $fusion = $x[0].';'.$y[0].';'.$x[1].';'.$y[1];
	        
	    #
	    # Extract the type of rearrangement between fusion partners
	    #
	    my($rearrangement_type,$distance) = extract_rearrangement_class($x[0],$y[0],$ORDER,$AORDER); # Import HUGO id
#	    print "$fusion\t$rearrangement_type\t$distance\t";
	    
	    #
	    # WORK IN PROGRESS (are they single best or multi-mapping; what is the exon coordinates; 
	    # 
	    
	    # Tally fusion
	    $FINAL{$fusion}[0]++;
	    $FINAL{$fusion}[1] = $$UCSC_INFO{$x[0]}[0];
	    $FINAL{$fusion}[2] = $$UCSC_INFO{$y[0]}[0];
	    
#	    print "$FINAL{$fusion}[0]\t$FINAL{$fusion}[1]\t$FINAL{$fusion}[2]\n";
	    
	    # Store positions within transcript 
	    push @{$FPpositions{$fusion}}, $x[3];
	    push @{$TPpositions{$fusion}}, $y[3];
	    
	    # Exon starts
	    push @{$FPexonstarts{$fusion}}, $x[6];
	    push @{$TPexonstarts{$fusion}}, $y[6];
	    
	    # Exon ends
	    push @{$FPexonends{$fusion}}, $x[7];
	    push @{$TPexonends{$fusion}}, $y[7];

	    # Store exons
	    push @{$FPexons{$fusion}}, $x[5];
	    push @{$TPexons{$fusion}}, $y[5];

	    # Store mate pair ids
	    push @{$MP{$fusion}}, $mp;
	    next;
	}
    }    
}


#
# Print finalized list of chimeras
#
my @chimeras = keys %FINAL;
foreach my $chimera (@chimeras){
#    print "$chimera\n";
    my($FPhugo,$TPhugo,$FPucsc,$TPucsc) = split(/\;/,$chimera); 
    #
    # Print those with greater than 2 reads
    #
    if($FINAL{$chimera}[0] >= '1'){
	
        #
	# Determine the proximity within the sample to determine if they fit predicted breakpoint
	#
	my $mps = @{$MP{$chimera}};
	#print "$chimera\t$mps\n";
	my($FPdistance,$FPscore) = determine_virtual_breakpoint(\@{$FPpositions{$chimera}});
	my($TPdistance,$FPscore) = determine_virtual_breakpoint(\@{$TPpositions{$chimera}});


	my($FPexon_range,$FPstr,$relstart,$FPtranscript_length);
####	if($FPucsc eq 'uc010vof.1' && $TPucsc eq 'uc002fir.1'){
	    my($FPexon_range) = extract_exons(\@{$FPexons{$chimera}}); 
	    my($FPES,$FPEE) = split(/\-/,$FPexon_range);
#	    print "$FPucsc\t$FPES\t$FPEE\n";
	    if($FPEE eq ''){
		($relstart,$FPtranscript_length,$FPstr) = extract_rel_exon('5p',$FPES,$FPucsc,$EXONS,$TPucsc);
#		if($FPucsc eq 'uc010gor.2' && $TPucsc eq 'uc010gnx.2'){
#                    print "HUH\t$relstart\t$FPtranscript_length\t$FPstr\n";
#                }
	    }else{
		($relstart,$FPtranscript_length,$FPstr) = extract_rel_exon('5p',$FPEE,$FPucsc,$EXONS,$TPucsc);
#		if($FPucsc eq 'uc010gor.2' && $TPucsc eq 'uc010gnx.2'){
#		    print "WHAT\t$relstart\t$FPtranscript_length\t$FPstr\n";
#		}
	    }
#	}

	
	my($TPexon_range) = extract_exons(\@{$TPexons{$chimera}}); 
	my($TPES,$TPEE) = split(/\-/,$TPexon_range); 
	my($relend,$TPtranscript_length,$TPstr) = extract_rel_exon('3p',$TPES,$TPucsc,$EXONS,$FPucsc);
	
#	if($FPucsc eq 'uc010gor.2' && $TPucsc eq 'uc010gnx.2'){
#	    print "FLIP\t$FPdistance\t$TPdistance\t$relstart\t$FPtranscript_length\t$FPstr\t";
#	    print "$relend\t$TPtranscript_length\t$TPstr\n";
#	}



	if($FPdistance ne 'Fail' && $TPdistance ne 'Fail'){
	    # if( ($FPdistance + $TPdistance <= '500') || $chimera eq 'TMPRSS2-ERG'){
	    print OUT "$chimera\t$FINAL{$chimera}[0]\t$FINAL{$chimera}[1]\t$FINAL{$chimera}[2]\t$FPdistance\t$TPdistance\t$FPexon_range\t$TPexon_range\t$FPstr\t$TPstr\t";
#	    print OUTBEDPE "$FPucsc\t$relstart\t$FPtranscript_length\t$TPucsc\t$relend\t$TPtranscript_length\t$FPhugo-$TPhugo\t$FINAL{$chimera}[0]\t$FPstr\t$TPstr\t";
#	    if($FPucsc eq 'uc010gor.2' && $TPucsc eq 'uc010gnx.2'){
		print "$FPucsc\t$relstart\t$FPtranscript_length\t$TPucsc\t$relend\t$TPtranscript_length\t$FPhugo-$TPhugo\t";
		print "$FINAL{$chimera}[0]\t$FPstr\t$TPstr\t";
#	    }

	    
	    #
	    # Print out 5' partner reads
	    #
	    foreach my $mp (@{$MP{$chimera}}){
		print "$SEQS5P{$mp}[0];";
	    }
	    print "\t";
	    foreach my $mp (@{$MP{$chimera}}){
		print "$SEQS3P{$mp}[0];";
	    }
	    print "\t";

	    foreach my $mp (@{$MP{$chimera}}){
		print "$mp,";
		#print OUT "$chimera\t$mp\n";
		print OUT "\t$mp";
		#my(@w)=split(/\;/,$SEQS{$mp}[0]);
		print READS ">".$chimera."_".$mp."_Mate1\n".$SEQS5P{$mp}[0]."\n>".$chimera."_".$mp."_Mate2\n".$SEQS3P{$mp}[0]."\n";
	    }
	    print OUT "\n";
	    print "\n";
	    # }
####	}
	}
    }
}
close(OUTBEDPE);

############################################
#                                          #
#                                          #
#             SUBROUTINES                  #
#                                          #
#                                          #
############################################
sub determine_position{
    my($partner,$exons,$ch,$s1,$e1,$exon_lengths,$exon_starts,$ucscstr,$mate_start,$mate_end,$ucsc) = @_;
    
    my $genomic_anchor = $s1;
    my $transcript_length ='0';
    #my $exon = '1';
    
    for(my $x = '0'; $x < $exons; $x++){
       	my $exon_start = $s1 + $$exon_starts[$x];
	my $exon_end = $s1 + $$exon_starts[$x] + $$exon_lengths[$x];
	my $exon_length = $exon_end - $exon_start;
       	
	#if($ucsc eq 'uc002yxb.1'){
	#    print "$partner\t$ucsc\t$ucscstr\t$x\t$mate_start\t$mate_end\tEXON\t$exon_start\t$exon_end\n";
	#}

        # If mate is 5p partner and UCSC transcript if positively oriented
	if($partner eq '5p' && $ucscstr eq '+'){
	    if($mate_start >= $exon_start && $mate_start <= $exon_end){
		#print "$x\t$ch\t$exon_start\t$exon_end\t$mate_start\n";
		my $exon_portion = $mate_start - $exon_start;
		my $mate_transcript_position = $exon_portion + $transcript_length;
		my $current_exon = $x + 1;
		return($mate_transcript_position,$current_exon,$exon_start,$exon_end);
	    }
	}

	# If mate is 5p partner and UCSC transcript if negatively oriented
	if($partner eq '5p' && $ucscstr eq '-'){
	    if($mate_end >= $exon_start && $mate_end <= $exon_end){
		my $exon_portion = $mate_end - $exon_start;
		my $mate_transcript_position = $exon_portion + $transcript_length;
		#if($ucsc eq 'uc002yzj.1'){
		#    print "TMPRSS2\t$mate_transcript_position\n";
		#}
		my $current_exon = $exons - $x;
		return($mate_transcript_position,$current_exon,$exon_start,$exon_end);
	    }
	}

	# If mate is 3p partner and UCSC transcript if negatively oriented
	if($partner eq '3p' && $ucscstr eq '-'){
	    if($mate_start >= $exon_start && $mate_start <= $exon_end){
		#print "$x\t$ch\t$exon_start\t$exon_end\t$mate_start\n";
		my $exon_portion = $mate_start - $exon_start;
		my $mate_transcript_position = $exon_portion + $transcript_length;
		#if($ucsc eq 'uc002yxb.1'){
		#    print "YOWSWER\t$mate_transcript_position\n";
		#}
		my $current_exon = $exons - $x;
                return($mate_transcript_position,$current_exon,$exon_start,$exon_end);
	    }
	}

	# If mate is 3p partner and UCSC transcript if positively oriented
	if($partner eq '3p' && $ucscstr eq '+'){
	    if($mate_end >= $exon_start && $mate_end <= $exon_end){
		my $exon_portion = $mate_end - $exon_start;
		my $mate_transcript_position = $exon_portion + $transcript_length;
		my $current_exon = $x + 1;
                return($mate_transcript_position,$current_exon,$exon_start,$exon_end);
	    }
	}
	$transcript_length += $exon_length;
    }    
}

sub extract_rel_exon{
    my($fusion_partner,$exon_of_interest,$ucsc,$exons_href,$other_gene) = @_;
    
    my $transcript_length = '0';
    my $point_of_interest = '0';

    my $chr = $$exons_href{$ucsc}[0]; #$v[1], $v[2], $v[7], $v[8], $v[9], $v[10]; # Chr, Str, exons, starts, ends, hugo 
    my $str = $$exons_href{$ucsc}[1];
    my $exons = $$exons_href{$ucsc}[2];
    my $starts = $$exons_href{$ucsc}[3]; my @starts = split(/\,/,$starts);# my $extra = pop(@starts);
    my $ends = $$exons_href{$ucsc}[4]; my @ends = split(/\,/,$ends); # my $extra = pop(@ends);
    my $hugo = $$exons_href{$ucsc}[5];
    
    my $exon_count = '1';
    if($str eq '+'){
	for(my $x = '0'; $x < $exons; $x++){
	    my $exon_start = $starts[$x];
	    my $exon_end = $ends[$x];
	    my $exon_length = $ends[$x] - $starts[$x];
#	    if($ucsc eq 'uc002yzj.2' && $fusion_partner eq '5p'){
		my $rel_exon_start = $transcript_length;
		my $rel_exon_end = $transcript_length + $exon_length;
		if($exon_count eq $exon_of_interest && $fusion_partner eq '5p'){
		    $point_of_interest = $rel_exon_end;
		}elsif($exon_count eq $exon_of_interest && $fusion_partner eq '3p'){
		    $point_of_interest = $rel_exon_start;
		}
	#	print "Exon\t$fusion_partner\t$x\t$exon_count\t$exon_of_interest\t$exon_start\t$exon_end\t$exon_length\t$ucsc\t$other_gene\t$str\t$rel_exon_start\t$rel_exon_end\t$point_of_interest\t$transcript_length\n";
		$transcript_length += $exon_length;
		$exon_count++;
	#    }
	}
    }elsif($str eq '-'){
	for(my $x = $exons-1; $x >= 0; $x--){
	    my $exon_start = $starts[$x];
	    my $exon_end = $ends[$x];
	    my $exon_length = $ends[$x] - $starts[$x];
#	    if($ucsc eq 'uc010gor.2' && $fusion_partner eq '5p' && $other_gene eq 'uc010gnx.2'){
	    my $rel_exon_start = $transcript_length;
	    my $rel_exon_end = $transcript_length + $exon_length;
	    if($exon_count eq $exon_of_interest && $fusion_partner eq '5p'){
		$point_of_interest = $rel_exon_end;
	    }elsif($exon_count eq $exon_of_interest && $fusion_partner eq '3p'){
		$point_of_interest = $rel_exon_start;
	    }
	#    print "Exon\t$fusion_partner\t$x\t$exon_count\t$exon_of_interest\t$exon_start\t$exon_end\t$exon_length\t$ucsc\t$other_gene\t$str\t$rel_exon_start\t$rel_exon_end\t$point_of_interest\t$transcript_length\n";
	    $transcript_length += $exon_length;
	    $exon_count++;
#	    }
	}
    }
#    if($ucsc eq 'uc010gor.2' && $fusion_partner eq '5p' && $other_gene eq 'uc010gnx.2'){
#	print "$ucsc\t$chr\t$point_of_interest\t$transcript_length\n";
#    }
    return($point_of_interest,$transcript_length,$str);
}


sub extract_exons{
    my($exons) = @_;

    my %tmp;
    # Remove null values
    foreach my $exon (@$exons){
        if($exon ne ''){
            $tmp{$exon}[0]++;
        }
    }
    
    my @exons = keys %tmp;
    my @sorted = sort { $a <=> $b } @exons;
    my $first = shift(@sorted);
    my $last = pop(@sorted);
    my $exon_range = $first.'-'.$last;
    return($exon_range);
}

sub determine_virtual_breakpoint{
    my($positions) = @_;
    
    my %tmp;
    # Remove null values
    foreach my $pos (@$positions){
	if($pos ne ''){
	    $tmp{$pos}[0]++;
	}
    }

    # Uniq supporting mappings
    my @uniq_pos = keys %tmp;
    my $score = @uniq_pos;
    
    # Sort array
    if($score >= '2'){
	my @sorted = sort { $a <=> $b } @uniq_pos;
	#foreach my $sorted (@sorted){
	#    print "$sorted ";
	#}print "\n";
	my $first = shift(@sorted);
	my $last = pop(@sorted);
	my $distance = $last - $first;
	#print "WOOO\t$first\t$last\t$distance\n";
	return($distance,$score);
    }else{
	return('Fail',$score);
    }
}

sub extract_rearrangement_class{
    my($FPgene,$TPgene,$ORDER,$AORDER) = @_;
    
    my $rearrangement_type;
    my $distance = 'N/A';
    my($type,$orientation,$overlapping);
    my($fchr,$fstr,$fst,$fend) = split(/\;/,$$UCSC_INFO{$FPgene}[0]);
    my($tchr,$tstr,$tst,$tend) = split(/\;/,$$UCSC_INFO{$TPgene}[0]);
    
    # Orientation
    if($fstr ne $tstr){
	$orientation = 'Opposite';
    }else{
	$orientation = 'Same';
    }

    $overlapping = '0';
    
    # Check for inter-chromosomal rearrangements
    if($tchr ne $fchr){
	return ("Inter-chromosomal",$distance);
    }else{
	
	# Overlapping
	if( ($tst >= $fst && $tst <= $fend) || ($tend >= $fst && $tend <= $fend) ){
	    $overlapping++;
	    #$rearrangement_type = 'Overlapping';
	    if($orientation eq 'Opposite' && $fst < $tst && $fstr eq '+'){
		return ('Overlapping_Converging',$distance);
	    }elsif($orientation eq 'Opposite' && $fst > $tst && $fstr eq '+'){
		return ('Overlapping_Diverging',$distance);
	    }elsif($orientation eq 'Same' && $fst < $tst && $fstr eq '+'){
		return ('Overlapping_Same',$distance);
	    }elsif($orientation eq 'Same' && $fst > $tst && $fstr eq '+'){
		return ('Overlapping_Complex',$distance);
	    }elsif($orientation eq 'Opposite' && $fst > $tst && $fstr eq '-'){
		return ('Overlapping_Converging',$distance);
	    }elsif($orientation eq 'Opposite' && $fst < $tst && $fstr eq '-'){
		return ('Overlapping_Diverging',$distance);
	    }elsif($orientation eq 'Same' && $fst > $tst && $fstr eq '-'){
		return ('Overlapping_Same',$distance);
	    }elsif($orientation eq 'Same' && $fst < $tst && $fstr eq '-'){
		return ('Overlapping_Complex',$distance);
	    }
	}

	# Sort out complex intra-rearrangements and read-throughs	
	
	if($fst <= $tst && $fstr eq '+'){
	    $distance = $tst - $fst;
	    #$FPgene = $fusion;
	}elsif($tst <= $fst && $fstr eq '+'){
	    $distance = $fst - $tst;
	    #$FPgene = $fusion;
	}elsif($fst >= $tst && $fstr eq '-'){
	    #$FPgene = $fusion;
	    $distance = $tst - $fst;
	}elsif($tst >= $fst && $fstr eq '-'){
	    #$FPgene = $fusion;
	    $distance = $fst - $tst;
	}
	
	my $tdistance = sprintf ("%.2f",$distance / 1000);
	$distance = $tdistance;
	
	if($orientation eq 'Same'){
	    my $genes_between_same_strand = check_for_readthrough($fchr,$fstr,$FPgene,$TPgene,$ORDER);
	    if($genes_between_same_strand <= '2'){
		return ('Read-through',$distance);
	    }else{
		return ('Intra-chromosomal_read',$distance);
	    }
	}elsif($orientation eq 'Opposite'){
	    my $genes_between = check_for_neighbors($fchr,$FPgene,$TPgene,$AORDER);
	    
	    #if($FPgene eq 'ARFGEF2'){
	#	print "$FPgene\t$TPgene\t$genes_between\t$distance\t$tst\t$tend\t$tstr\t$fst\t$fend\t$fstr\t$orientation\n";
	#    }
	    
	    if($genes_between <= '2'){
		if($fst < $tst && $fst eq '+'){
		    return ('Adjacent_Converging',$distance);
		}elsif($fst > $tst && $fstr eq '+'){
		    return ('Adjacent_Diverging',$distance);
		}elsif($fst < $tst && $fstr eq '+'){
		    return ('Adjacent_Same',$distance);
		}elsif($fst > $tst && $fstr eq '+'){
		    return ('Adjacent_Complex',$distance);
		}elsif($fst > $tst && $fstr eq '-'){
		    return ('Adjacent_Converging',$distance);
		}elsif($fst < $tst && $fstr eq '-'){
		    return ('Adjacent_Diverging',$distance);
		}elsif($fst > $tst && $fstr eq '-'){
		    return ('Adjacent_Same',$distance);
		}elsif($fst < $tst && $fstr eq '-'){
		    return ('Adjacent_Complex',$distance);
		}
	    }else{
		if($fst < $tst && $fstr eq '+'){
		    return ('Intra-chromosomal_complex',$distance);
		}elsif($fst > $tst && $fstr eq '+'){
		    return ('Intra-chromosomal_complex',$distance);
		}elsif($fst < $tst && $fstr eq '+'){
		    return ('Intra-chromosomal_complex',$distance);
		}elsif($fst > $tst && $fstr eq '+'){
		    return ('Intra-chromosomal_complex',$distance);
		}elsif($fst > $tst && $fstr eq '-'){
		   return ('Intra-chromosomal_complex',$distance);
		}elsif($fst < $tst && $fstr eq '-'){
		    return ('Intra-chromosomal_complex',$distance);
		}elsif($fst > $tst && $fstr eq '-'){
		    return ('Intra-chromosomal_complex',$distance);
		}elsif($fst < $tst && $fstr eq '-'){
		    return ('Intra-chromosomal_complex',$distance);
		}
	    }
	}
    }
    return ('Undetermined',$distance);
}
    
sub check_for_readthrough{
    my($chr,$strand,$FPgene,$TPgene,$ORDER) = @_;
    
    my $neighbors = '0';
    my $chr_strand = $chr.'_'.$strand;
    
    # Need to make sure the genes in between the fusion genes are in fact independent
    my %TMP;
    my $store_genes = '0';
    foreach my $gene (@{$$ORDER{$chr_strand}}){
	if($strand eq '+'){
	    if($store_genes eq '1' && $gene eq $TPgene){
		last;
	    }
	    
	    if($store_genes eq '1' && $gene ne $TPgene){
		if($gene ne $TPgene && $gene ne $FPgene){
		    $TMP{$gene}[0]++;
		}
		next;
	    }
	    
	    if($gene eq $FPgene){
		$store_genes = '1';
		next;
	    }
	}else{
	    if($store_genes eq '1' && $gene eq $FPgene){
		last;
	    }
	    
	    if($store_genes eq '1' && $gene ne $FPgene){
		if($gene ne $TPgene && $gene ne $FPgene){
		    $TMP{$gene}[0]++;
		}
		next;
	    }
	    
	    if($gene eq $TPgene){
		$store_genes = '1';
		next;
	    }
	}	
    }	
    
    # Process the TMP hash
    my @genes = keys %TMP;
    my $genes_in_between_same_orientation = @genes;

    return($genes_in_between_same_orientation);
}

sub check_for_neighbors{
    my($chr,$FPgene,$TPgene,$AORDER) = @_;
    
    my $neighbors = '0';
    
    # Need to make sure the genes in between the fusion genes are in fact independent
    my %TMP;
    my $store_gene1 = '0';
    my $store_gene2 = '0';
    foreach my $gene (@{$$AORDER{$chr}}){
	if($store_gene1 eq '1' && $gene eq $TPgene){
	    last;
	}
	
	if($store_gene2 eq '1' && $gene eq $FPgene){
	    last;
	}
	
	if($store_gene1 eq '1' && $gene ne $TPgene){
	    if($gene ne $TPgene && $gene ne $FPgene){
		$TMP{$gene}[0]++;
	    }
	    next;
	}

	if($store_gene2 eq '1' && $gene ne $FPgene){
	    if($gene ne $TPgene && $gene ne $FPgene){
		$TMP{$gene}[0]++;
	    }
	    next;
	}
	
	if($gene eq $FPgene){
	    $store_gene1 = '1';
	    next;
	}
	
	if($gene eq $TPgene){
	    $store_gene2 = '1';
	    next;
	}
    }	
    
    # Process the TMP hash
    my @genes = keys %TMP;
    my $genes_in_between = @genes;

    return($genes_in_between);
}
    
sub get_ucsc{
    my($ucsc_table) = @_;

    my(%UCSC,%UCSC2GENE,%UCSCID,%ORDER,%AORDER,%EXONS);
    open(UCSC, "<$ucsc_table");
    while(<UCSC>){
        chomp;
		
		if(/^\#/){
	    	next;
		}else{
	   	 #my($UCSC,$chr,$str,$tst,$tend,$cst,$cend,$exons,$gene)=split(/\t/);
	   	 my(@v) = split(/\t/);
	   	 my $id = $v[1].';'.$v[2].';'.$v[3].';'.$v[4];
	   	 push @{$UCSC{$v[10]}}, $id; 
	   	 $UCSC2GENE{$v[0]}[0] = $v[10];
	   	 $UCSCID{$v[0]}[0] = $id;
	   	 my $index = $v[1].'_'.$v[2];
	   	 push @{$ORDER{$index}}, $v[10];
	   	 push @{$AORDER{$v[1]}}, $v[10];
	   	 push @{$EXONS{$v[0]}}, $v[1], $v[2], $v[7], $v[8], $v[9], $v[10]; # Chr, Str, exons, starts, ends, hugo 
		}
    }
    close(UCSC);
    
    return(\%UCSC,\%UCSC2GENE,\%UCSCID,\%ORDER,\%AORDER,\%EXONS);
}
###################################################################################
#
# Christopher Maher
# chrmaher@med.mich.edu
# Last modified : 03/27/2008
#
###################################################################################
