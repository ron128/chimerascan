'''
Created on Feb 25, 2012

@author: mkiyer
'''
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

def build_tid_gene_map(bamfh, genefile, rname_prefix=None):
    rname_tid_map = dict((rname,tid) for tid,rname in enumerate(bamfh.references))
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

def build_transcript_cluster_map(line_iter, rname_prefix=None):
    # setup cluster trees
    chrom_strand_cluster_trees = \
        collections.defaultdict(lambda: {"+": ClusterTree(0,1),
                                         "-": ClusterTree(0,1)})
    transcripts = []
    index_cluster_map = {}
    for transcript in GeneFeature.parse(line_iter):
        # insert exons into cluster tree
        cluster_tree = chrom_strand_cluster_trees[transcript.chrom][transcript.strand]
        i = len(transcripts)
        for start,end in transcript.exons:
            cluster_tree.insert(start, end, i)
        # each transcript is initially in a cluster by itself
        index_cluster_map[i] = set([i])
        transcripts.append(transcript)
    # extract gene clusters
    for strand_cluster_trees in chrom_strand_cluster_trees.itervalues():
        for cluster_tree in strand_cluster_trees.itervalues():
            for start, end, indexes in cluster_tree.getregions():
                # make new cluster by aggregating all existing
                # clusters with new indexes
                newclust = set(indexes)
                for i in indexes:
                    newclust.update(index_cluster_map[i])
                # map every transcript to the new cluster
                for i in newclust:
                    index_cluster_map[i] = newclust
    # enumerate all clusters
    rname_prefix = '' if rname_prefix is None else rname_prefix
    transcript_cluster_map = {}
    for cluster_id, clust in enumerate(index_cluster_map.values()):
        for i in clust:
            transcript = transcripts[i]
            transcript_cluster_map[rname_prefix + transcript.tx_name] = cluster_id
    return transcript_cluster_map

def build_transcript_tid_cluster_map(bamfh, line_iter, rname_prefix=None):
    # make the standard cluster map
    transcript_cluster_map = build_transcript_cluster_map(line_iter, rname_prefix)
    # map reference name to tid
    transcript_tid_map = {}
    rname_prefix = '' if rname_prefix is None else rname_prefix
    for tid,rname in enumerate(bamfh.references):
        if rname.startswith(rname_prefix):
            transcript_tid_map[rname] = tid
    # remake the cluster map
    tid_cluster_map = {}
    for rname, cluster_id in transcript_cluster_map.iteritems():
        if rname not in transcript_tid_map:
            continue
        tid = transcript_tid_map[rname]
        tid_cluster_map[tid] = cluster_id
    return tid_cluster_map

def build_transcript_genome_map(line_iter, rname_prefix=None):
    # create arrays to map genes in bed file to genome 
    rname_prefix = '' if rname_prefix is None else rname_prefix
    transcript_genome_map = {}    
    for g in GeneFeature.parse(line_iter):
        rname = rname_prefix + g.tx_name
        strand = 1 if g.strand == '-' else 0 
        exon_vectors = [(start, end) for start, end in g.exons]
        if strand:
            exon_vectors.reverse()
        if rname in transcript_genome_map:
            logging.error("Duplicate references %s found in bed file" % (rname))
        transcript_genome_map[rname] = (g.chrom, strand, exon_vectors)
    return transcript_genome_map

def build_transcript_tid_genome_map(bamfh, line_iter, rname_prefix=None):
    # make the standard map
    transcript_genome_map = build_transcript_genome_map(line_iter, rname_prefix)
    # map reference name to tid
    rname_prefix = '' if rname_prefix is None else rname_prefix
    transcript_tid_map = {}
    for tid,rname in enumerate(bamfh.references):
        if rname.startswith(rname_prefix):
            transcript_tid_map[rname] = tid
    # remap using tid as key
    tid_genome_map = {}
    for rname, coords in transcript_genome_map.iteritems():
        if rname not in transcript_tid_map:
            continue
        tid = transcript_tid_map[rname]
        tid_genome_map[tid] = coords
    return tid_genome_map

def transcript_to_genome_pos(rname, pos, transcript_genome_map):    
    '''
    translate gene 'rname' position 'gene_pos' to genomic
    coordinates.  returns a 3-tuple with (chrom, strand, pos)
    '''
    chrom, strand, intervals = transcript_genome_map[rname]
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