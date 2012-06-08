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

def build_tid_transcript_map(bamfh, feature_iter):
    rname_tid_map = dict((rname,tid) for tid,rname in enumerate(bamfh.references))
    tid_tx_map = {}
    # build gene and genome data structures for fast lookup
    for f in feature_iter:
        tid = rname_tid_map[str(f.tx_id)]
        tid_tx_map[tid] = f
    return tid_tx_map

def build_transcript_map(feature_iter):
    tx_map = {}
    # build gene and genome data structures for fast lookup
    for f in feature_iter:
        tx_map[str(f.tx_id)] = f
    return tx_map

def build_genome_transcript_trees(feature_iter):
    genome_tx_trees = collections.defaultdict(lambda: IntervalTree())    
    # build gene and genome data structures for fast lookup
    for g in feature_iter:
        # add gene to interval tree
        interval = Interval(g.tx_start, g.tx_end, strand=g.strand, value=g)
        genome_tx_trees[g.chrom].insert_interval(interval)
    return genome_tx_trees

def cluster_transcripts(feature_iter):
    # setup cluster trees
    chrom_strand_cluster_trees = \
        collections.defaultdict(lambda: {"+": ClusterTree(0,1),
                                         "-": ClusterTree(0,1)})
    transcripts = []
    index_cluster_map = {}
    for feature in feature_iter:
        # insert exons into cluster tree
        cluster_tree = chrom_strand_cluster_trees[feature.chrom][feature.strand]
        i = len(transcripts)
        for start,end in feature.exons:
            cluster_tree.insert(start, end, i)
        # each transcript is initially in a cluster by itself
        index_cluster_map[i] = set([i])
        transcripts.append(feature)
    # extract gene clusters
    for strand_cluster_trees in chrom_strand_cluster_trees.itervalues():
        for cluster_tree in strand_cluster_trees.itervalues():
            for start, end, indexes in cluster_tree.getregions():
                # make new cluster by aggregating all existing
                # clusters with new indexes
                newclust = set(indexes)
                for i in indexes:
                    newclust.update(index_cluster_map[i])
                newclust = frozenset(newclust)
                # map every transcript to the new cluster
                for i in newclust:
                    index_cluster_map[i] = newclust
    # consolidate list of clusters
    clusters = set()
    for clust in index_cluster_map.itervalues():
        if clust not in clusters:
            clusters.add(clust)
    for clust in sorted(clusters):
        yield tuple(transcripts[i] for i in clust)

def build_transcript_cluster_map(feature_iter):
    transcript_cluster_map = {}
    for cluster_id, transcripts in enumerate(cluster_transcripts(feature_iter)):
        for t in transcripts:        
            transcript_cluster_map[str(t.tx_id)] = cluster_id
    return transcript_cluster_map

def build_transcript_tid_cluster_map(bamfh, feature_iter):
    # make the standard cluster map
    transcript_cluster_map = build_transcript_cluster_map(feature_iter)
    # map reference name to tid
    transcript_tid_map = {}
    for tid,rname in enumerate(bamfh.references):
        transcript_tid_map[rname] = tid
    # remake the cluster map
    tid_cluster_map = {}
    for rname, cluster_id in transcript_cluster_map.iteritems():
        if rname not in transcript_tid_map:
            continue
        tid = transcript_tid_map[rname]
        tid_cluster_map[tid] = cluster_id
    return tid_cluster_map

def build_transcript_genome_map(feature_iter):
    # create arrays to map genes in bed file to genome 
    transcript_genome_map = {}    
    for f in feature_iter:
        strand = 1 if f.strand == '-' else 0 
        exon_vectors = [(start, end) for start, end in f.exons]
        if strand:
            exon_vectors.reverse()
        tx_id = str(f.tx_id)            
        if tx_id in transcript_genome_map:
            logging.error("Duplicate references %s found in bed file" % (tx_id))
        transcript_genome_map[tx_id] = (f.chrom, strand, exon_vectors)
    return transcript_genome_map

def build_tid_transcript_genome_map(bamfh, feature_iter):
    # make the standard map
    transcript_genome_map = build_transcript_genome_map(feature_iter)
    # map reference name to tid
    transcript_tid_map = {}
    for tid,rname in enumerate(bamfh.references):
        transcript_tid_map[rname] = tid
    # remap using tid as key
    tid_genome_map = {}
    for rname, coords in transcript_genome_map.iteritems():
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