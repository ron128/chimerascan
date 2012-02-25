'''
Created on Jan 31, 2011

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
import collections

from chimerascan.bx.cluster import ClusterTree
from chimerascan.bx.intersection import Interval, IntervalTree
# local imports
from feature import GeneFeature

def get_rname_tid_map(bamfh):
    rname_tid_map = {}
    for tid,ref in enumerate(bamfh.references):
        rname_tid_map[ref] = tid
    return rname_tid_map

def build_tid_gene_map(bamfh, genefile, rname_prefix=None):
    rname_tid_map = get_rname_tid_map(bamfh)
    rname_prefix = '' if rname_prefix is None else rname_prefix
    tid_tx_map = {}
    # build gene and genome data structures for fast lookup
    for g in GeneFeature.parse(open(genefile)):
        # only use genes that are references in the sam file
        rname = rname_prefix + g.tx_name
        if rname not in rname_tid_map:
            continue
        tid = rname_tid_map[rname]
        tid_tx_map[tid] = g
    return tid_tx_map

def build_tx_name_gene_map(genefile, rname_prefix=None):
    rname_prefix = '' if rname_prefix is None else rname_prefix
    tx_map = {}
    # build gene and genome data structures for fast lookup
    for g in GeneFeature.parse(open(genefile)):
        tx_map[rname_prefix + g.tx_name] = g
    return tx_map

def build_genome_tx_trees(genefile):
    genome_tx_trees = collections.defaultdict(lambda: IntervalTree())    
    # build gene and genome data structures for fast lookup
    for g in GeneFeature.parse(open(genefile)):
        # add gene to interval tree
        interval = Interval(g.tx_start, g.tx_end, strand=g.strand, value=g)
        genome_tx_trees[g.chrom].insert_interval(interval)
    return genome_tx_trees

def build_rname_cluster_map(line_iter, rname_prefix=None):
    rname_prefix = '' if rname_prefix is None else rname_prefix
    cluster_trees = collections.defaultdict(lambda: ClusterTree(0,1))
    genes = []    
    for g in GeneFeature.parse(line_iter):
        rname = rname_prefix + g.tx_name
        # insert into cluster tree        
        cluster_trees[g.chrom].insert(g.tx_start, g.tx_end, len(genes)) 
        genes.append(g)
    # extract gene clusters
    tx_cluster_map = {}
    current_cluster_id = 0
    for chrom, tree in cluster_trees.iteritems():
        for start, end, indexes in tree.getregions():
            # group overlapping transcripts on same strand together            
            strand_tx_dict = collections.defaultdict(lambda: set())
            for index in indexes:
                g = genes[index]
                rname = rname_prefix + g.tx_name              
                strand_tx_dict[g.strand].add(rname)
            # build a map between transcript tids and all the overlapping
            # transcripts on the same strand
            for strand, rnames in strand_tx_dict.iteritems():
                for rname in rnames:                    
                    tx_cluster_map[rname] = current_cluster_id
                current_cluster_id += 1
    return tx_cluster_map

def build_rname_genome_map(line_iter, rname_prefix=None):
    # create arrays to map genes in bed file to genome 
    rname_prefix = '' if rname_prefix is None else rname_prefix
    gene_genome_map = {}    
    for g in GeneFeature.parse(line_iter):
        rname = rname_prefix + g.tx_name
        strand = 1 if g.strand == '-' else 0 
        exon_vectors = [(start, end) for start, end in g.exons]
        if strand:
            exon_vectors.reverse()
        if rname in gene_genome_map:
            logging.error("Duplicate references %s found in bed file" % (rname))
        gene_genome_map[rname] = (g.chrom, strand, exon_vectors)
    return gene_genome_map

def transcript_to_genome_pos(rname, pos, gene_genome_map):    
    '''
    translate gene 'rname' position 'gene_pos' to genomic
    coordinates.  returns a 3-tuple with (chrom, strand, pos)
    '''
    chrom, strand, intervals = gene_genome_map[rname]
    offset = 0
    for start, end, in intervals:
        exon_size = end - start
        if pos < offset + exon_size:            
            if strand:
                return chrom, strand, start + exon_size - (pos - offset) - 1
            else:
                return chrom, strand, start + (pos - offset)
        #print start, end, offset, pos
        offset += exon_size
    return None