'''
Created on Jul 21, 2011

@author: mkiyer

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

import pysam

from chimerascan.lib import config
from chimerascan.lib.feature import TranscriptFeature
from chimerascan.lib.chimera import DiscordantRead, Chimera, frags_to_encomp_string
from chimerascan.lib.transcriptome import build_transcript_map
from chimerascan.lib.fragment_size_distribution import InsertSizeDistribution
from chimerascan.lib.seq import calc_homology

def parse_discordant_bedpe_by_transcript_pair(fh):
    prev_tx5p, prev_tx3p = None,None
    frags = []
    for line in fh:
        fields = line.strip().split('\t')        
        tx5p = fields[0]
        tx3p = fields[3]
        dr5p = DiscordantRead.from_list(fields[10].split("|"))
        dr3p = DiscordantRead.from_list(fields[11].split("|"))
        if (tx5p, tx3p) != (prev_tx5p, prev_tx3p):
            if len(frags) > 0:
                yield prev_tx5p, prev_tx3p, frags
                frags = []
            prev_tx5p, prev_tx3p = tx5p, tx3p
        frags.append((dr5p, dr3p))
    if len(frags) > 0:
        yield tx5p, tx3p, frags        

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
    return best_isize_prob, best_breakpoints

def extract_breakpoint_sequence(tx_id_5p, tx_end_5p, 
                                tx_id_3p, tx_start_3p, 
                                ref_fa, max_read_length,
                                homology_mismatches):
    tx_start_5p = max(0, tx_end_5p - max_read_length + 1)
    tx_end_3p = tx_start_3p + max_read_length - 1
    # fetch sequence
    seq5p = ref_fa.fetch(tx_id_5p, tx_start_5p, tx_end_5p).upper()
    seq3p = ref_fa.fetch(tx_id_3p, tx_start_3p, tx_end_3p).upper()
    # pad sequence if too short
    if len(seq5p) < (max_read_length - 1):
        logging.warning("Could not extract sequence of length >%d from "
                        "5' partner at %s:%d-%d, only retrieved "
                        "sequence of length %d" % 
                        (max_read_length-1, tx_id_5p, tx_start_5p, 
                         tx_end_5p, len(seq5p)))
        # pad sequence
        padding = (max_read_length - 1) - len(seq5p)
        seq5p = ("N" * padding) + seq5p
    if len(seq3p) < max_read_length - 1:
        logging.warning("Could not extract sequence of length >%d from "
                        "3' partner at %s:%d-%d, only retrieved "
                        "sequence of length %d" % 
                        (max_read_length-1, tx_id_3p, tx_start_3p, 
                         tx_end_3p, len(seq3p)))
        # pad sequence
        padding = (max_read_length - 1) - len(seq3p)
        seq3p = seq3p + ("N" * padding)
    # if 5' partner continues along its normal transcript
    # without fusing, get the sequence that would result
    homolog_end_5p = tx_end_5p + max_read_length - 1
    homolog_seq_5p = ref_fa.fetch(tx_id_5p, tx_end_5p, homolog_end_5p).upper()
    # if 3' partner were to continue in the 5' direction,
    # grab the sequence that would be produced
    homolog_start_3p = max(0, tx_start_3p - max_read_length + 1)
    homolog_seq_3p = ref_fa.fetch(tx_id_3p, homolog_start_3p, tx_start_3p).upper()
    # count number of bases in common between downstream 5' sequence
    # and the sequence of the 3' partner in the chimera
    homology_right = calc_homology(homolog_seq_5p, seq3p, 
                                   homology_mismatches)
    # count number of bases in common between upstream 3' sequence
    # and the sequence of the 5' partner in the chimera
    homology_left = calc_homology(homolog_seq_3p[::-1], seq5p[::-1],
                                  homology_mismatches)
    return seq5p, seq3p, homology_left, homology_right

def nominate_chimeras(index_dir, isize_dist_file, input_file, output_file, 
                      trim_bp, max_read_length, homology_mismatches):
    # read insert size distribution
    isize_dist = InsertSizeDistribution.from_file(open(isize_dist_file))
    # build a lookup table to get genomic intervals from transcripts
    logging.debug("Reading transcript information")
    transcript_feature_file = os.path.join(index_dir, config.TRANSCRIPT_FEATURE_FILE)    
    transcripts = list(TranscriptFeature.parse(open(transcript_feature_file)))
    tx_id_map = build_transcript_map(transcripts)
    # open the reference sequence fasta file
    ref_fasta_file = os.path.join(index_dir, config.TRANSCRIPTOME_FASTA_FILE)
    ref_fa = pysam.Fastafile(ref_fasta_file)
    # keep track of mapping from breakpoint sequence to breakpoint id
    # this requires storing all breakpoint sequences in memory which is
    # potentially expensive.  TODO: investigate whether this should be
    # moved to a separate sort-update-sort procedure
    breakpoint_seq_name_map = {}
    breakpoint_num = 1
    # group discordant read pairs by gene
    logging.debug("Parsing discordant reads")
    chimera_num = 1
    outfh = open(output_file, "w")    
    for tx_id_5p, tx_id_3p, frags in parse_discordant_bedpe_by_transcript_pair(open(input_file)):
        # get gene information
        tx5p = tx_id_map[tx_id_5p]
        tx3p = tx_id_map[tx_id_3p]
        # bin fragments into putative breakpoints
        breakpoint_dict = collections.defaultdict(lambda: [])
        for dr5p,dr3p in frags:
            # given the insert size find the highest probability 
            # exon junction breakpoint between the two transcripts
            isize_prob, breakpoints = \
                choose_best_breakpoints(dr5p, dr3p, tx5p, tx3p, 
                                        trim_bp, isize_dist)
            for breakpoint in breakpoints:
                breakpoint_dict[breakpoint].append((dr5p, dr3p))        
        # iterate through breakpoints and build chimera candidates
        for breakpoint,frags in breakpoint_dict.iteritems():          
            exon_num_5p, tx_end_5p, exon_num_3p, tx_start_3p = breakpoint
            breakpoint_seq_5p, breakpoint_seq_3p, homology_left, homology_right = \
                extract_breakpoint_sequence(tx_id_5p, tx_end_5p,
                                            tx_id_3p, tx_start_3p,
                                            ref_fa, max_read_length,
                                            homology_mismatches)                
            tx3p_length = sum((end - start) for start,end in tx3p.exons)
            # get unique breakpoint id based on sequence
            breakpoint_seq = breakpoint_seq_5p + breakpoint_seq_3p
            if breakpoint_seq in breakpoint_seq_name_map:
                breakpoint_name = breakpoint_seq_name_map[breakpoint_seq]
            else:
                breakpoint_name = "B%07d" % (breakpoint_num)
                breakpoint_seq_name_map[breakpoint_seq] = breakpoint_name
                breakpoint_num += 1
            # write gene, breakpoint, and raw reads to a file and follow the
            # BEDPE format
            gene_names_5p = ",".join(sorted(set(["_".join(x.split()) for x in tx5p.gene_names])))
            gene_names_3p = ",".join(sorted(set(["_".join(x.split()) for x in tx3p.gene_names])))
            fields = [tx5p.tx_id, 0, tx_end_5p,  # chrom1, start1, end1
                      tx3p.tx_id, tx_start_3p, tx3p_length, # chrom2, start2, end2
                      "C%07d" % (chimera_num), # name
                      1.0, # pvalue
                      tx5p.strand, tx3p.strand, # strand1, strand2
                      gene_names_5p, gene_names_3p, # gene names
                      # exon interval information
                      '%d-%d' % (0, exon_num_5p),
                      '%d-%d' % (exon_num_3p, len(tx3p.exons)),
                      # breakpoint information
                      breakpoint_name, 
                      breakpoint_seq_5p, breakpoint_seq_3p, 
                      homology_left, homology_right, 
                      # fragments
                      frags_to_encomp_string(frags),
                      # spanning reads
                      None]
            print >>outfh, '\t'.join(map(str, fields))
            chimera_num += 1
    outfh.close()
    ref_fa.close()
    return config.JOB_SUCCESS
    

def main():
    from optparse import OptionParser
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = OptionParser("usage: %prog [options] <index> <isize_dist.txt> "
                          "<discordant_reads.srt.bedpe> <chimeras.txt>")
    parser.add_option("--trim", dest="trim", type="int", 
                      default=config.EXON_JUNCTION_TRIM_BP,
                      help="apply trimming when choosing exon boundaries to "
                           "to consider possible breakpoints")
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
    input_file = args[2]
    output_file = args[3]
    return nominate_chimeras(index_dir, isize_dist_file, 
                             input_file, output_file, 
                             options.trim,
                             options.max_read_length,
                             options.homology_mismatches)


if __name__ == '__main__':
    sys.exit(main())
