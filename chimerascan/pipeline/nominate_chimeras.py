'''
Created on May 26, 2011

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
import collections

from chimerascan import pysam
from chimerascan.bx.cluster import ClusterTree

from chimerascan.lib import config
from chimerascan.lib.base import FR_UNSTRANDED, FR_FIRSTSTRAND, \
    FR_SECONDSTRAND, NO_STRAND, POS_STRAND, NEG_STRAND
from chimerascan.lib.sam import get_strand, get_genomic_intervals, parse_pe_reads
from chimerascan.lib.transcriptome import build_exon_interval_trees, get_predicted_strand
from chimerascan.lib.chimera import DiscordantRead

def imax2(a,b):
    return a if a >= b else b

def get_read_strand_info(r, chrom, intervals, library_type, exon_trees):
    strand = get_strand(r)
    if library_type == FR_UNSTRANDED:
        align_strand = NEG_STRAND if r.is_reverse else POS_STRAND
        if strand == NO_STRAND:
            # try to infer strand from overlapping transcripts
            strand = get_predicted_strand(chrom, intervals, exon_trees)
            if strand == NO_STRAND:
                # still don't know strand?
                sense = None
        else:
            sense = (align_strand == strand)        
    else:
        # strand is known because this is a strand-specific protocol
        if library_type == FR_FIRSTSTRAND:
            # dUTP method: read1 is antisense and read2 is sense
            sense = r.is_read2
        elif library_type == FR_SECONDSTRAND:
            sense = r.is_read1
    return strand, sense

def resolve_unknown_sense_info(pe_dreads):
    """
    for reads with unknown sense orientation, try to use the paired 
    mate read's orientation (if known) to predict the correct 
    orientation
    """
    r1_dirs = set(d.sense for d in pe_dreads[0])
    r2_dirs = set(d.sense for d in pe_dreads[1])
    r1_has_unknown = False
    r2_has_unknown = False
    if None in r1_dirs:
        r1_has_unknown = True
        r1_dirs.remove(None)
    if None in r2_dirs:
        r2_has_unknown = True
        r2_dirs.remove(None)
    # if there is only one orientation (5' or 3')
    # then can use it to resolve the other orientation 
    if r1_has_unknown and len(r2_dirs) == 1:
        r2_sense = list(r2_dirs)[0]
        # TODO: to support 'ff' libraries, need to adjust
        # this code to check library type before assuming
        # read is on opposite strand of mate
        r1_sense = (not r2_sense)        
        for d in pe_dreads[0]:
            if d.sense is None:
                d.sense = r1_sense
    if r2_has_unknown and len(r1_dirs) == 1:
        r1_sense = list(r1_dirs)[0]
        r2_sense = (not r1_sense)        
        for d in pe_dreads[1]:
            if d.sense is None:
                d.sense = r2_sense

def get_discordant_read_info(bamfh, r, inner_dist, library_type, exon_trees):
    d = DiscordantRead()
    d.qname = r.qname
    d.hit_index=r.opt('HI')
    d.chrom = bamfh.getrname(r.rname)
    d.readnum = 0 if r.is_read1 else 1
    intervals = get_genomic_intervals(r)
    d.strand, d.sense = get_read_strand_info(r, d.chrom, intervals, 
                                                library_type, 
                                                exon_trees)
    # extend right from rightmost interval of read
    d.right_interval = (intervals[-1][0],
                        intervals[-1][1] + inner_dist)
    # extend left from leftmost interval of read
    d.left_interval = (imax2(0, intervals[0][0] - inner_dist),
                       intervals[0][1])
    return d


def pair_discordant_reads(dr1, dr2):
    if (dr1.sense is None) and (dr2.sense is None):
        # neither of the alignments have sense information
        # TODO: handle this
        return False, None
    elif dr1.sense is None:
        # check read2 for information and infer sense
        sense1 = not dr2.sense
        sense2 = dr2.sense
    elif dr2.sense is None:
        sense1 = dr1.sense
        sense2 = not dr1.sense
    elif dr1.sense == dr2.sense:
        # one read must be sense and the other must be
        # antisense
        return False, None
    if dr1.sense:
        interval5p = dr1.right_interval if (dr1.strand == POS_STRAND) else dr1.left_interval
        interval3p = dr2.right_interval if (dr2.strand == NEG_STRAND) else dr2.left_interval                 
        return True, dr1, dr2, interval5p, interval3p
    else:
        interval5p = dr2.right_interval if (dr1.strand == POS_STRAND) else dr2.left_interval
        interval3p = dr1.right_interval if (dr2.strand == NEG_STRAND) else dr1.left_interval            
        return True, dr2, dr1, interval5p, interval3p

def localize_breakpoints(bamfile, gene_file, inner_dist, library_type):
    # read transcript information and build interval trees
    exon_trees = build_exon_interval_trees(gene_file)
    # parse discordant reads
    bamfh = pysam.Samfile(bamfile, "rb")
    for pe_reads in parse_pe_reads(bamfh):
        # determine breakpoint intervals and orientation
        # of discordant reads
        pe_dreads = [[], []]
        for readnum, reads in enumerate(pe_reads):
            pe_reads[readnum] = [get_discordant_read_info(bamfh, r, inner_dist, 
                                                          library_type, exon_trees) 
                                 for r in reads]
        for dr1 in pe_dreads[0]:
            for dr2 in pe_dreads[1]:
                
        
        # 
        # if orientation of one of the reads is known, can use it to 
        # adjust orientation of the other read
        #resolve_unknown_sense_info(pe_dreads)

#    trees = [collections.defaultdict(lambda: ClusterTree(0, 1)),
#             collections.defaultdict(lambda: ClusterTree(0, 1))]            
#    # cluster breakpoint regions
#    for readnum, dinfos in enumerate(pe_dreads):
#        for d in dinfos:
#            dir = 0 if d.sense else 1
#            id = len(discordant_infos)
#            if d.strand == NO_STRAND:
#                if d.sense is None:
#                    dirs = (0, 1)
#                
#            if ((d.strand == POS_STRAND and d.sense) or
#                (d.strand == NEG_STRAND and (not d.sense))):                    
#                start, end = d.right_start, d.right_end
#                trees[dir][d.chrom].insert(d.right_start, d.right_end, id)
#            else:
                    
                
        
    bamfh.close()


def get_breakpoint_interval(chrom, intervals, strand, sense,
                            inner_dist):

    if ((strand == POS_STRAND and sense) or
        (strand == NEG_STRAND and (not sense))):
        # 5' partner is on "+" stranded gene, OR
        # 3' partner is on "-" strand 
        # so extend right from rightmost interval of read
        start = intervals[-1][0]
        end = intervals[-1][1] + inner_dist
    elif ((strand == POS_STRAND and (not sense)) or
          (strand == NEG_STRAND and sense)):
        # 3' partner is on "+" stranded gene, so extend left 
        # from leftmost interval of read
        start = imax2(0, intervals[0][0] - inner_dist)
        end = intervals[0][1]
    return (chrom, start, end)


def nominate_chimeras(bamfile, gene_file, max_fragment_length):
    bamfh = pysam.Samfile(bamfile, "rb")    
    bamfh.close()



def main():
    from optparse import OptionParser
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = OptionParser("usage: %prog [options] <index> <segment_length> <discordant_reads.bam> <chimeras.bedpe>")
    parser.add_option('--max-fragment-length', dest="max_fragment_length", 
                      type="int", default=1000)
    parser.add_option('--library-type', dest="library_type")
    options, args = parser.parse_args()
    index_dir = args[0]
    segment_length = int(args[1])
    bamfile = args[2]

    gene_file = os.path.join(index_dir, config.GENE_FEATURE_FILE)

    input_file = args[0]
    output_file = args[1]
    gene_file = os.path.join(options.index_dir, config.GENE_FEATURE_FILE)
    nominate_chimeras(open(input_file), 
                      open(output_file,'w'),
                      gene_file)


if __name__ == '__main__':
    main()
