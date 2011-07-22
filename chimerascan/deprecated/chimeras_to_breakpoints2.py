'''
Created on Jul 15, 2011

@author: mkiyer

chimerascan: chimeric transcript discovery using RNA-seq

Copyright (C) 2011 Matthew Iyer

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''
import logging
import os
import sys
import collections
import itertools
import operator

from chimerascan import pysam

from chimerascan.lib import config
from chimerascan.lib.chimera import Chimera, \
    DiscordantTags, DISCORDANT_TAG_NAME, \
    OrientationTags, ORIENTATION_TAG_NAME, \
    DiscordantRead, ChimeraTypes, ChimeraPartner
from chimerascan.lib.breakpoint import Breakpoint
from chimerascan.lib.gene_to_genome import build_tid_tx_map
from chimerascan.lib.fragment_size_distribution import InsertSizeDistribution
from chimerascan.lib.seq import calc_homology

def parse_pairs(bamfh):
    bam_iter = iter(bamfh)
    try:
        while True:
            r1 = bam_iter.next()
            r2 = bam_iter.next()
            yield r1,r2
    except StopIteration:
        pass
    
def parse_gene_discordant_reads(bamfh):
    """
    return tuples of (5',3') reads that both align to transcripts
    """
    for r1,r2 in parse_pairs(bamfh):
        # TODO:
        # for now we are only going to deal with gene-gene
        # chimeras and leave other chimeras for study at a 
        # later time
        dr1 = r1.opt(DISCORDANT_TAG_NAME)
        dr2 = r2.opt(DISCORDANT_TAG_NAME)
        if (dr1 != DiscordantTags.DISCORDANT_GENE or
            dr2 != DiscordantTags.DISCORDANT_GENE):            
            continue
        # organize key in 5' to 3' order
        or1 = r1.opt(ORIENTATION_TAG_NAME)
        or2 = r2.opt(ORIENTATION_TAG_NAME)
        assert or1 != or2
        if or1 == OrientationTags.FIVEPRIME:
            pair = (r1,r2)
        else:
            pair = (r2,r1)
        yield pair

def get_chimera_type(fiveprime_gene, threeprime_gene, gene_trees):
    """
    return tuple containing ChimeraType and distance 
    between 5' and 3' genes 
    """
    # get gene information
    chrom5p, start5p, end5p, strand1 = fiveprime_gene.chrom, fiveprime_gene.tx_start, fiveprime_gene.tx_end, fiveprime_gene.strand
    chrom3p, start3p, end3p, strand2 = threeprime_gene.chrom, threeprime_gene.tx_start, threeprime_gene.tx_end, threeprime_gene.strand
    # interchromosomal
    if chrom5p != chrom3p:
        return ChimeraTypes.INTERCHROMOSOMAL, None
    # orientation
    same_strand = strand1 == strand2
    # genes on same chromosome so check overlap
    is_overlapping = (start5p < end3p) and (start3p < end5p)            
    if is_overlapping:
        if not same_strand:
            if ((start5p <= start3p and strand1 == "+") or
                (start5p > start3p and strand1 == "-")):                    
                return (ChimeraTypes.OVERLAP_CONVERGE, 0)
            else:
                return (ChimeraTypes.OVERLAP_DIVERGE, 0)
        else:
            if ((start5p <= start3p and strand1 == "+") or
                (end5p >= end3p and strand1 == "-")):
                return (ChimeraTypes.OVERLAP_SAME, 0)
            else:
                return (ChimeraTypes.OVERLAP_COMPLEX, 0)
    # if code gets here then the genes are on the same chromosome but do not
    # overlap.  first calculate distance (minimum distance between genes)
    if start5p <= start3p:
        distance = start3p - end5p
        between_start,between_end = end5p,start3p
    else:
        distance = end3p - start5p
        between_start,between_end = end3p,start5p
    # check whether there are genes intervening between the
    # chimera candidates
    genes_between = []
    genes_between_same_strand = []
    for hit in gene_trees[chrom5p].find(between_start,
                                       between_end):
        if (hit.start > between_start and
            hit.end < between_end):             
            if hit.strand == strand1:
                genes_between_same_strand.append(hit)
            genes_between.append(hit)
            
    if same_strand:
        if len(genes_between_same_strand) == 0:
            return ChimeraTypes.READTHROUGH, distance
        else:
            return ChimeraTypes.INTRACHROMOSOMAL, distance
    else:
        # check for reads between neighboring genes    
        if len(genes_between) == 0:
            if ((start5p <= start3p and strand1 == "+") or
                (start5p > start3p and strand1 == "-")):                    
                return (ChimeraTypes.ADJ_CONVERGE, distance)
            elif ((start5p >= start3p and strand1 == "+") or
                  (start5p < start3p and strand1 == "-")):
                return (ChimeraTypes.ADJ_DIVERGE, distance)
            elif ((start5p <= start3p and strand1 == "+") or
                  (start5p > start3p and strand1 == "-")):
                return (ChimeraTypes.ADJ_SAME, distance)
            elif ((start5p >= start3p and strand1 == "+") or
                  (start5p < start3p and strand1 == '-')):
                return (ChimeraTypes.ADJ_COMPLEX, distance)
        else:
            return ChimeraTypes.INTRA_COMPLEX, distance    
    return ChimeraTypes.UNKNOWN, distance


def read_pairs_to_chimera(chimera_name, tid5p, tid3p, readpairs, 
                          tid_tx_map, genome_tx_trees, trim_bp):
    # get gene information
    tx5p = tid_tx_map[tid5p]
    tx3p = tid_tx_map[tid3p]
    # categorize chimera type
    chimera_type, distance = get_chimera_type(tx5p, tx3p, genome_tx_trees)
    # create chimera object
    c = Chimera()
    iter5p = itertools.imap(operator.itemgetter(0), readpairs)
    iter3p = itertools.imap(operator.itemgetter(1), readpairs)
    c.partner5p = ChimeraPartner.from_discordant_reads(iter5p, tx5p, trim_bp)
    c.partner3p = ChimeraPartner.from_discordant_reads(iter3p, tx3p, trim_bp)
    c.name = chimera_name
    c.chimera_type = chimera_type
    c.distance = distance
    # raw reads
    c.encomp_read_pairs = readpairs
    return c

def calc_isize_prob(isize, isize_dist):
    # find percentile of observing this insert size in the reads
    isize_per = isize_dist.percentile_at_isize(isize)
    # convert to a probability score (0.0-1.0)
    isize_prob = 1.0 - (2.0 * abs(50.0 - isize_per))/100.0    
    return isize_prob

def choose_best_breakpoints(r5p, r3p, tx5p, tx3p, trim_bp, isize_dist):
    best_breakpoints = set()
    best_isize_prob = None
    # iterate through 5' transcript exons    
    exon_iter_5p = reversed(tx5p.exons) if tx5p.strand == '-' else iter(tx5p.exons)
    tx_end_5p = 0
    for exon_num_5p,coords5p in enumerate(exon_iter_5p):
        genome_start_5p, genome_end_5p = coords5p        
        exon_size_5p = genome_end_5p - genome_start_5p
        tx_end_5p += exon_size_5p
        # fast forward on 5' gene to first exon beyond read        
        if tx_end_5p < (r5p.aend - trim_bp):
            continue        
        #print "tx end 5p", tx_end_5p, "exon_size_5p", exon_size_5p, "r5p.aend", r5p.aend, "trim_bp", trim_bp
        # now have a candidate insert size between between 5' read and
        # end of 5' exon
        isize5p = tx_end_5p - r5p.pos
        # iterate through 3' transcript
        exon_iter_3p = reversed(tx3p.exons) if tx3p.strand == '-' else iter(tx3p.exons)
        tx_start_3p = 0
        local_best_breakpoints = set()
        local_best_isize_prob = None
        for exon_num_3p,coords3p in enumerate(exon_iter_3p):
            genome_start_3p, genome_end_3p = coords3p
            #print "\t", coords3p 
            # stop after going past read on 3' transcript
            if tx_start_3p >= (r3p.pos + trim_bp):
                break
            # get another candidate insert size between start of 3'
            # exon and 3' read
            isize3p = r3p.aend - tx_start_3p
            #print "\t", isize5p, isize3p, tx_end_5p, tx_start_3p
            # compare the insert size against the known insert size
            # distribution
            isize_prob = calc_isize_prob(isize5p + isize3p, isize_dist)
            if ((local_best_isize_prob is None) or
                (isize_prob > local_best_isize_prob)):
                local_best_isize_prob = isize_prob
                local_best_breakpoints = set([(exon_num_5p, tx_end_5p, 
                                               exon_num_3p, tx_start_3p)])
            elif (isize_prob == local_best_isize_prob):
                local_best_breakpoints.add((exon_num_5p, tx_end_5p, 
                                            exon_num_3p, tx_start_3p))
            tx_start_3p += genome_end_3p - genome_start_3p
        # compare locally best insert size probability to global best
        if ((best_isize_prob is None) or 
            (local_best_isize_prob > best_isize_prob)):
            best_isize_prob = local_best_isize_prob
            best_breakpoints = local_best_breakpoints
        elif (local_best_isize_prob == best_isize_prob):
            # for ties we keep all possible breakpoints
            best_breakpoints.update(local_best_breakpoints)
    # TODO: remove debugging output
    #ends5p = [x[1] for x in best_breakpoints]
    #starts3p = [x[3] for x in best_breakpoints]
    #print ends5p, starts3p, "r1:%d-%d" % (r5p.pos, r5p.aend), "r2:%d-%d" % (r3p.pos, r3p.aend), best_isize_prob
    return best_isize_prob, best_breakpoints

def extract_breakpoint_sequence(tx_name_5p, tx_end_5p, 
                                tx_name_3p, tx_start_3p, 
                                ref_fa, max_read_length,
                                homology_mismatches):
    tx_start_5p = max(0, tx_end_5p - max_read_length + 1)
    tx_end_3p = tx_start_3p + max_read_length - 1
    # fetch sequence
    seq5p = ref_fa.fetch(tx_name_5p, tx_start_5p, tx_end_5p)
    seq3p = ref_fa.fetch(tx_name_3p, tx_start_3p, tx_end_3p)
    # pad sequence if too short
    if len(seq5p) < (max_read_length - 1):
        logging.warning("Could not extract sequence of length >%d from "
                        "5' partner at %s:%d-%d, only retrieved "
                        "sequence of length %d" % 
                        (max_read_length-1, tx_name_5p, tx_start_5p, 
                         tx_end_5p, len(seq5p)))
        # pad sequence
        padding = (max_read_length - 1) - len(seq5p)
        seq5p = ("N" * padding) + seq5p
    if len(seq3p) < max_read_length - 1:
        logging.warning("Could not extract sequence of length >%d from "
                        "3' partner at %s:%d-%d, only retrieved "
                        "sequence of length %d" % 
                        (max_read_length-1, tx_name_3p, tx_start_3p, 
                         tx_end_3p, len(seq3p)))
        # pad sequence
        padding = (max_read_length - 1) - len(seq3p)
        seq3p = seq3p + ("N" * padding)
    # if 5' partner continues along its normal transcript
    # without fusing, get the sequence that would result
    homolog_end_5p = tx_end_5p + max_read_length - 1
    homolog_seq_5p = ref_fa.fetch(tx_name_5p, tx_end_5p, homolog_end_5p)
    # if 3' partner were to continue in the 5' direction,
    # grab the sequence that would be produced
    homolog_start_3p = max(0, tx_start_3p - max_read_length + 1)
    homolog_seq_3p = ref_fa.fetch(tx_name_3p, homolog_start_3p, tx_start_3p)
    # count number of bases in common between downstream 5' sequence
    # and the sequence of the 3' partner in the chimera
    homology_right = calc_homology(homolog_seq_5p, seq3p, 
                                   homology_mismatches)
    # count number of bases in common between upstream 3' sequence
    # and the sequence of the 5' partner in the chimera
    homology_left = calc_homology(homolog_seq_3p[::-1], seq5p[::-1],
                                  homology_mismatches)
    return seq5p, seq3p, homology_left, homology_right

def discordant_reads_to_breakpoints(index_dir, isize_dist_file, 
                                    input_bam_file, output_file, 
                                    trim_bp, max_read_length,
                                    homology_mismatches):                      
    """
    homology_mismatches: number of mismatches to tolerate while computing
    homology between chimeric breakpoint sequence and "wildtype" sequence
    
    trim_bp: when selecting the best matching exon for each read, we
    account for spurious overlap into adjacent exons by trimming the
    read by 'trim_bp'
    """   
    # read insert size distribution
    isize_dist = InsertSizeDistribution.from_file(open(isize_dist_file))
    # open BAM alignment file
    bamfh = pysam.Samfile(input_bam_file, "rb")
    # build a lookup table to get genomic intervals from transcripts
    logging.debug("Reading gene information")
    gene_file = os.path.join(index_dir, config.GENE_FEATURE_FILE)
    tid_tx_map = build_tid_tx_map(bamfh, gene_file,
                                  rname_prefix=config.GENE_REF_PREFIX)
    # open the reference sequence fasta file
    ref_fasta_file = os.path.join(index_dir, config.ALIGN_INDEX + ".fa")
    ref_fa = pysam.Fastafile(ref_fasta_file)
    # iterate through read pairs
    outfh = open(output_file, "w")
    logging.debug("Parsing discordant reads")
    for r5p,r3p in parse_gene_discordant_reads(bamfh):
        # store pertinent read information in lightweight structure called
        # DiscordantRead object. this departs from SAM format into a 
        # custom read format
        dr5p = DiscordantRead.from_read(r5p)
        dr3p = DiscordantRead.from_read(r3p)
        # get gene information
        tx5p = tid_tx_map[r5p.rname]
        tx3p = tid_tx_map[r3p.rname]
        # given the insert size find the highest probability 
        # exon junction breakpoint between the two transcripts
        isize_prob, breakpoints = \
            choose_best_breakpoints(r5p, r3p, tx5p, tx3p, 
                                    trim_bp, isize_dist)
        # extract the sequence of the breakpoint along with the
        # number of homologous bases at the breakpoint between 
        # chimera and wildtype genes
        for breakpoint in breakpoints:
            exon_num_5p, tx_end_5p, exon_num_3p, tx_start_3p = breakpoint
            breakpoint_seq_5p, breakpoint_seq_3p, homology_left, homology_right = \
                extract_breakpoint_sequence(config.GENE_REF_PREFIX + tx5p.tx_name, tx_end_5p,
                                            config.GENE_REF_PREFIX + tx3p.tx_name, tx_start_3p,
                                            ref_fa, max_read_length,
                                            homology_mismatches)
            # write breakpoint information for each read to a file
            fields = [tx5p.tx_name, 0, tx_end_5p,
                      tx3p.tx_name, tx_start_3p, tx3p.tx_end,
                      r5p.rname,  # name
                      isize_prob, # score
                      tx5p.strand, tx3p.strand, # strand 1, strand 2
                      # user defined fields
                      exon_num_5p, exon_num_3p,
                      breakpoint_seq_5p, breakpoint_seq_3p, 
                      homology_left, homology_right] 
            fields.append('|'.join(map(str, dr5p.to_list())))
            fields.append('|'.join(map(str, dr3p.to_list())))  
            print >>outfh, '\t'.join(map(str, fields))        
    # cleanup
    ref_fa.close()
    outfh.close()
    bamfh.close()
    return config.JOB_SUCCESS

def main():
    from optparse import OptionParser
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = OptionParser("usage: %prog [options] <index> <isizedist.txt> "
                          "<discordant_reads.bedpe> <chimeras.txt>")
    parser.add_option("--trim", dest="trim", type="int", 
                      default=config.EXON_JUNCTION_TRIM_BP, metavar="N",
                      help="Trim ends of reads by N bp when determining "
                      "the start/end exon of chimeras [default=%default]")
    parser.add_option("--max-read-length", dest="max_read_length", type="int",
                      default=100, metavar="N",
                      help="Reads in the BAM file are guaranteed to have "
                      "length less than N [default=%default]")
    parser.add_option("--homology-mismatches", type="int", 
                      dest="homology_mismatches", 
                      default=config.BREAKPOINT_HOMOLOGY_MISMATCHES,
                      help="Number of mismatches to tolerate when computing "
                      "homology between gene and its chimeric partner "
                      "[default=%default]")
    options, args = parser.parse_args()
    index_dir = args[0]
    isize_dist_file = args[1]
    input_bam_file = args[2]
    output_file = args[3]
    return nominate_chimeras(index_dir, isize_dist_file, 
                             input_bam_file, output_file, 
                             trim_bp=options.trim,
                             max_read_length=options.max_read_length,
                             homology_mismatches=options.homology_mismatches)


if __name__ == '__main__':
    sys.exit(main())
