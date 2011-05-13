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

from feature import GeneFeature
from sam import get_genomic_intervals
from chimerascan.bx.intersection import Interval, IntervalTree

def build_exon_interval_trees(genefile):
    exon_trees = collections.defaultdict(lambda: IntervalTree)
    exon_intervals = {}
    # build gene and genome data structures for fast lookup
    for g in GeneFeature.parse(open(genefile)):
        for i,e in enumerate(g.exons):
            k = (g.chrom, e[0], e[1])
            if k not in exon_intervals:
                # add exon to tree
                txlist = []
                exon_intervals[k] = txlist
                exon_trees[g.chrom].insert_interval(e[0], e[1], 
                                                    strand=g.strand, 
                                                    value=exon_intervals[k])
            else:
                txlist = exon_intervals[k]
            # add transcript isoform
            txlist.append((g,i))
    return exon_intervals, exon_trees

def get_transcripts_at_interval(chrom, start, end, exon_intervals, exon_trees):
    hits = []
    for hit in exon_trees[chrom].find(start, end):
        txlist = hit.value
        # check for compatibility with overlapping genes
        for g,exon_num in txlist:
            pass


def get_transcript_coords(read, exon_intervals, exon_trees):
    intervals = get_genomic_intervals(read)
    # get all transcripts compatible with first interval
    
    
    pass



