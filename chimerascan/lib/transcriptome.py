'''
Created on May 2, 2011

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
import collections
import logging

from base import cmp_strand, NO_STRAND
from feature import GeneFeature
from sam import get_genomic_intervals, get_strand
from chimerascan.bx.intersection import Interval, IntervalTree

def build_gene_interval_trees(genefile):
    trees = collections.defaultdict(lambda: IntervalTree())
    intervals = {}
    # build gene interval trees for fast lookup by genomic position
    for g in GeneFeature.parse(open(genefile)):
        k = (g.chrom, g.tx_start, g.tx_end)
        if k not in intervals:
            # add gene to tree
            txlist = []
            intervals[k] = txlist
            interval = Interval(g.tx_start, g.tx_end, strand=g.strand, value=txlist)
            trees[g.chrom].insert_interval(interval)
        else:
            txlist = intervals[k]
        # add isoform to value (list of isoforms that share start/end)
        txlist.append(g)
    return trees

def get_genes_at_interval(chrom, start, end, strand, trees):
    features = []
    for hit in trees[chrom].find(start, end):
        txlist = hit.value
        # check for compatibility with overlapping genes
        # TODO: debug strand-specific libraries such that
        # antisense reads are not counted as part of overlapping
        # sense genes. this should work but needs testing
        features.extend(g for g in txlist if cmp_strand(g.strand, strand))
    return features

def get_overlapping_genes(bamfh, read, gene_trees):
    intervals = get_genomic_intervals(read)    
    # get read strand if it exists    
    strand = get_strand(read)
    # get all genes compatible with intervals
    hits = []
    for interval in intervals:
        hits.append(get_genes_at_interval(bamfh.getrname(read.rname),
                                          start=interval[0],
                                          end=interval[1],
                                          strand=strand,
                                          trees=gene_trees))
    return hits

def build_exon_interval_trees(genefile):
    exon_trees = collections.defaultdict(lambda: IntervalTree())
    exon_intervals = {}
    # build gene and genome data structures for fast lookup
    for g in GeneFeature.parse(open(genefile)):
        for i,e in enumerate(g.exons):
            k = (g.chrom, e[0], e[1])
            if k not in exon_intervals:
                # add exon to tree
                txlist = []
                exon_intervals[k] = txlist
                interval = Interval(e[0], e[1], strand=g.strand, value=txlist)
                exon_trees[g.chrom].insert_interval(interval)
            else:
                txlist = exon_intervals[k]
            # add transcript isoform
            txlist.append((g,i))
    return exon_trees


def get_transcripts_at_interval(chrom, start, end, strand, exon_trees):
    features = []
    for hit in exon_trees[chrom].find(start, end):
        txlist = hit.value
        for g,exon_num in txlist:
            features.extend(g for g in txlist if cmp_strand(g.strand, strand))
    return features

def get_overlapping_transcripts(bamfh, intervals, chrom, strand, exon_trees):
    hits = []
    for interval in intervals:
        hits.append(get_transcripts_at_interval(chrom=chrom,
                                                start=interval[0],
                                                end=interval[1],
                                                strand=strand,
                                                trees=exon_trees))
    return hits

def get_predicted_strand(bamfh, read, exon_trees):
    """predict strand of read by looking at overlapping genes"""
    intervals = get_genomic_intervals(read)
    hitlists = get_overlapping_transcripts(bamfh, 
                                           intervals, 
                                           bamfh.getrname(read.rname),
                                           NO_STRAND,
                                           exon_trees)
    strands = None
    for hits in hitlists:
        if strands is None:            
            # initialize with strand of transcripts at 
            # first interval
            strands = set(g.strand for g in hits)
        else:
            # intersect to find most compatible strands 
            strands.intersection_update(g.strand for g in hits)
    if (strands is None) or (len(strands) != 1):
        return NO_STRAND
    return strands.pop()

#class TranscriptAlignment(object):
#    def __init__(self, g, exon_num, chrom, start, end, e_start, e_end, 
#                 e_start_overhang, e_end_overhang, 
#                 tx_start, tx_end):
#        self.g = g
#        self.exon_num = exon_num
#        self.chrom = chrom
#        self.start = start
#        self.end = end
#        self.e_start = e_start
#        self.e_end = e_end
#        self.e_start_overhang = e_start_overhang
#        self.e_end_overhang = e_end_overhang
#        self.tx_start = tx_start
#        self.tx_end = tx_end
#
#    def __repr__(self):
#        return ("<%s(g='%s', exon_num='%d', chrom='%s', start='%d', "
#                "end='%d', e_start='%d', e_end='%d', "
#                "e_start_overhang='%d', e_end_overhang='%d', "
#                "tx_start='%d', tx_end='%d')>" %
#                (self.__class__.__name__, self.g, self.exon_num, 
#                 self.chrom, self.start, self.end, self.e_start, self.e_end, 
#                 self.e_start_overhang, self.e_end_overhang,
#                 self.tx_start, self.tx_end))
#
#def get_transcripts_at_interval(chrom, start, end, strand, exon_trees):
#    txsegs = []
#    for hit in exon_trees[chrom].find(start, end):
#        txlist = hit.value
#        # check for compatibility with overlapping genes
#        for g,exon_num in txlist:
#            g_exon_start, g_exon_end = g.exons[exon_num]            
#            #print 'ALN', g.gene_name, g.tx_name, 'enum', exon_num, g.strand, strand, g_exon_start, g_exon_end, start, end            
#            # strand must be compatible
#            if not cmp_strand(g.strand, strand):
#                continue
#            # get exon coordinates
#            g_exon_start, g_exon_end = g.exons[exon_num]            
#            e_start = max(start, g_exon_start) - g_exon_start
#            e_end = min(end, g_exon_end) - g_exon_start
#            e_start_overhang = max(0, g_exon_start - start)
#            e_end_overhang = max(0, end - g_exon_end)
#            # get transcript coordinates
#            tx_start = sum((end-start) for start,end in g.exons[:exon_num])
#            tx_start = tx_start + e_start
#            tx_end = tx_start + (end - start)
#            # convert to negative strand if necessary
#            if g.strand == NEG_STRAND:               
#                tx_length = sum((end - start) for start,end in g.exons)
#                tx_end, tx_start = (tx_length - tx_start, tx_length - tx_end)
#                exon_num = len(g.exons) - exon_num
#            txsegs.append(TranscriptAlignment(g=g, 
#                                              exon_num=exon_num,
#                                              chrom=chrom,
#                                              start=start,
#                                              end=end,                                               
#                                              e_start=e_start, 
#                                              e_end=e_end,
#                                              e_start_overhang=e_start_overhang,
#                                              e_end_overhang=e_end_overhang,
#                                              tx_start=tx_start,
#                                              tx_end=tx_end))
#    return txsegs
#
#
#def get_transcript_coords(bamfh, read, exon_intervals, exon_trees):
#    intervals = get_genomic_intervals(read)    
#    # get read strand if it exists    
#    strand = get_strand(read)
#    # get all transcripts compatible with first interval
#    txhits = []
#    for interval in intervals:
#        #print 'READ', read.is_read1, 'INTERVAL', interval
#        txhits.append(get_transcripts_at_interval(bamfh.getrname(read.rname),
#                                                  start=interval[0],
#                                                  end=interval[1],
#                                                  strand=strand,
#                                                  exon_intervals=exon_intervals,
#                                                  exon_trees=exon_trees))
#    return txhits
