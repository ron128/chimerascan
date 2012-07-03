'''
Created on Jun 20, 2012

@author: mkiyer
'''
import logging
import argparse
import sys
import os
import collections
import shelve

import pysam

from chimerascan.bx.cluster import ClusterTree
from chimerascan.lib import config
from chimerascan.lib.sam import get_aligned_intervals
from chimerascan.lib.chimera import ORIENTATION_TAG, ORIENTATION_5P, \
    ORIENTATION_3P, DISCORDANT_CLUSTER_TAG, DiscordantCluster, \
    discordant_cluster_to_string

def window_overlap(a, b):
    if a[0] != b[0]:
        return False
    return (a[1] <= b[2]) and (b[1] <= a[2])

def cluster_loci(bam_iter):
    try:
        # initialize window
        window = [bam_iter.next()]
        window_range = (window[0].tid, window[0].pos, window[0].aend)
        # separate into loci
        for r in bam_iter:
            # check if next read is outside current window
            interval = (r.tid, r.pos, r.aend)
            if not window_overlap(interval, window_range):
                # yield current window
                yield window
                # reset window
                window = [r]
                window_range = (r.tid, r.pos, r.aend)
            else:
                # add read to window
                window.append(r)
                window_range = (r.tid,
                                min(window_range[1], r.pos),
                                max(window_range[2], r.aend))
    except StopIteration:
        pass
    # yield last window
    if len(window) > 0:
        yield window

def get_concordant_frags(bamfh, rname, start, end, strand, orientation):
    # TODO: remove assert
    assert strand in ("+", "-")
    qnames = set()
    read_iter = bamfh.fetch(rname, start, end)
    # concordant fragments that extend past the cluster boundary are
    # evidence of non-chimeric transcripts 
    if (((strand == "+") and (orientation == ORIENTATION_5P)) or
        ((strand == "-") and (orientation == ORIENTATION_3P))):
        for r in read_iter:
            if (r.aend > end) or (r.pnext >= end):
                qnames.add(r.qname)
    else:
        for r in read_iter:
            if (r.pos < start) or (r.pnext < start):
                qnames.add(r.qname)
    return qnames

def get_unpaired_frags(bamfh, rname, start, end, strand, orientation):
    qnames = set()
    for r in bamfh.fetch(rname, start, end):
        # get genomic strand (+ or -) and orientation (5' or 3')
        rstrand = r.opt('XS')
        rorientation = r.opt(ORIENTATION_TAG)
        if (rstrand == strand) and (rorientation == orientation):
            qnames.add(r.qname)
    return qnames

def create_cluster(rname, start, end, cluster_id, strand, orientation, 
                   reads, unpaired_bamfh, concordant_bamfh):
    cluster_tree = ClusterTree(0,1)
    qnames = []
    for i,r in enumerate(reads):
        # add cluster tag to every read
        tagdict = collections.OrderedDict(r.tags)
        tagdict[DISCORDANT_CLUSTER_TAG] = cluster_id
        r.tags = tagdict.items()
        # keep read names
        qnames.append(r.qname)
        # cluster "exonic" intervals
        intervals = get_aligned_intervals(r)
        for istart,iend in intervals:
            cluster_tree.insert(istart, iend, i)
    # determine exon intervals
    exons = []
    for start, end, indexes in cluster_tree.getregions():
        exons.append((start,end))
    # count unpaired fragments where mapped mate aligns within cluster
    unpaired_qnames = get_unpaired_frags(unpaired_bamfh, 
                                         rname, start, end, 
                                         strand, orientation)
    # count wild-type non-chimeric fragments spanning cluster
    concordant_qnames = get_concordant_frags(concordant_bamfh, 
                                             rname, start, end, 
                                             strand, orientation)
    # make cluster object
    cluster = DiscordantCluster(rname=rname,
                                start=start,
                                end=end,
                                cluster_id=cluster_id,
                                strand=strand,
                                orientation=orientation,
                                exons=exons,
                                qnames=qnames,
                                unpaired_qnames=unpaired_qnames,                               
                                concordant_frags=len(concordant_qnames))
    return cluster

def add_reads_to_clusters(reads, next_cluster_id, discordant_bamfh, 
                          unpaired_bamfh, concordant_bamfh):
    # insert reads into clusters
    cluster_trees = {("+", ORIENTATION_5P): collections.defaultdict(lambda: ClusterTree(0,1)),
                     ("+", ORIENTATION_3P): collections.defaultdict(lambda: ClusterTree(0,1)),
                     ("-", ORIENTATION_5P): collections.defaultdict(lambda: ClusterTree(0,1)),
                     ("-", ORIENTATION_3P): collections.defaultdict(lambda: ClusterTree(0,1))}
    for i,r in enumerate(reads):
        # get genomic strand (+ or -) and orientation (5' or 3')
        strand = r.opt('XS')
        orientation = r.opt(ORIENTATION_TAG)
        # get appropriate cluster tree
        cluster_tree = cluster_trees[(strand, orientation)][r.tid]
        # insert read qname into tree
        cluster_tree.insert(r.pos, r.aend, i)
    # get read clusters
    clusters = []
    for strand_orientation, tid_cluster_trees in cluster_trees.iteritems():
        strand, orientation = strand_orientation
        for tid, cluster_tree in tid_cluster_trees.iteritems():
            rname = discordant_bamfh.getrname(tid)
            for start, end, indexes in cluster_tree.getregions():
                cluster_id = next_cluster_id
                cluster_reads = [reads[i] for i in indexes]
                cluster = create_cluster(rname, start, end, cluster_id, 
                                         strand, orientation, 
                                         cluster_reads, unpaired_bamfh, 
                                         concordant_bamfh)
                clusters.append(cluster)
                next_cluster_id += 1
    return clusters, next_cluster_id

def cluster_discordant_reads(discordant_bam_file, 
                             unpaired_bam_file,
                             concordant_bam_file, 
                             output_bam_file, 
                             cluster_file,
                             cluster_shelve_file):
    #
    # iterate through sorted discordant read alignments and form clusters
    # of overlapping alignments
    #
    logging.debug("Annotating discordant clusters")
    discordant_bamfh = pysam.Samfile(discordant_bam_file, "rb")
    unpaired_bamfh = pysam.Samfile(unpaired_bam_file, 'rb')
    concordant_bamfh = pysam.Samfile(concordant_bam_file, "rb")
    outbamfh = pysam.Samfile(output_bam_file, "wb", template=discordant_bamfh)
    outfh = open(cluster_file, "w")
    db = shelve.open(cluster_shelve_file)
    next_cluster_id = 0
    for locus_reads in cluster_loci(iter(discordant_bamfh)):
        locus_clusters, next_cluster_id = \
            add_reads_to_clusters(locus_reads, next_cluster_id, 
                                  discordant_bamfh, unpaired_bamfh, 
                                  concordant_bamfh)
        for cluster in locus_clusters:
            # write to shelve database
            db[str(cluster.cluster_id)] = cluster
            # write as tab-delimited text
            print >>outfh, discordant_cluster_to_string(cluster)
        # reads that now have cluster id tag so rewrite
        for r in locus_reads:
            outbamfh.write(r)
    db.close()
    outfh.close()
    outbamfh.close()
    concordant_bamfh.close()
    unpaired_bamfh.close()
    discordant_bamfh.close()
    #
    # index the newly annotated discordant bam file 
    #
    logging.debug("Indexing newly annotated discordant BAM file")
    pysam.index(output_bam_file)
    return config.JOB_SUCCESS
    
def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("discordant_bam_file") 
    parser.add_argument("unpaired_bam_file") 
    parser.add_argument("concordant_bam_file") 
    parser.add_argument("output_bam_file") 
    parser.add_argument("cluster_file")
    parser.add_argument("cluster_shelve_file")
    args = parser.parse_args()
    return cluster_discordant_reads(args.discordant_bam_file,
                                    args.unpaired_bam_file,
                                    args.concordant_bam_file, 
                                    args.output_bam_file, 
                                    args.cluster_file,
                                    args.cluster_shelve_file)

if __name__ == '__main__':
    sys.exit(main())