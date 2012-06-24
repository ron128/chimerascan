'''
Created on Jun 20, 2012

@author: mkiyer
'''
import logging
import argparse
import sys
import os
import collections

import pysam

from chimerascan.bx.cluster import ClusterTree

from chimerascan.lib import config
from chimerascan.lib.sam import parse_pe_reads
from chimerascan.lib.chimera import ORIENTATION_TAG, ORIENTATION_5P, \
    ORIENTATION_3P, DISCORDANT_CLUSTER_TAG, DiscordantCluster

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

def get_concordant_frags(bamfh, tid, start, end, strand, orientation):
    # TODO: remove assert
    assert strand in ("+", "-")
    qnames = set()
    rname = bamfh.getrname(tid)
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

def add_reads_to_clusters(reads, next_cluster_id, concordant_bamfh):
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
            for start, end, indexes in cluster_tree.getregions():
                cluster_id = next_cluster_id
                qnames = []
                for i in indexes:
                    r = reads[i]
                    qnames.append(r.qname)
                    # add cluster tag to every read
                    tagdict = collections.OrderedDict(r.tags)
                    tagdict[DISCORDANT_CLUSTER_TAG] = cluster_id
                    r.tags = tagdict.items()
                # count wild-type non-chimeric fragments spanning cluster
                concordant_qnames = get_concordant_frags(concordant_bamfh, 
                                                         tid, start, end, 
                                                         strand, orientation)
                cluster = DiscordantCluster(tid=tid,
                                            start=start,
                                            end=end,
                                            cluster_id=cluster_id,
                                            strand=strand,
                                            orientation=orientation,
                                            qnames=qnames,
                                            concordant_frags=len(concordant_qnames))
                clusters.append(cluster)
                next_cluster_id += 1
    return clusters, next_cluster_id

def cluster_discordant_reads(discordant_bam_file, 
                             concordant_bam_file, 
                             output_bam_file, 
                             cluster_file,
                             cluster_pair_file,
                             tmp_dir):
    if tmp_dir is None:
        tmp_dir = os.getcwd()
    #
    # iterate through sorted discordant read alignments and form clusters
    # of overlapping alignments
    #
    logging.debug("Annotating discordant clusters")
    discordant_bamfh = pysam.Samfile(discordant_bam_file, "rb")
    concordant_bamfh = pysam.Samfile(concordant_bam_file, "rb")
    outbamfh = pysam.Samfile(output_bam_file, "wb", template=discordant_bamfh)
    outfh = open(cluster_file, "w")
    next_cluster_id = 0
    for locus_reads in cluster_loci(iter(discordant_bamfh)):
        locus_clusters, next_cluster_id = \
            add_reads_to_clusters(locus_reads, next_cluster_id, 
                                  concordant_bamfh)
        for cluster in locus_clusters:
            fields = [discordant_bamfh.getrname(cluster.tid), 
                      cluster.start, cluster.end, cluster.cluster_id, 
                      cluster.strand, cluster.orientation, 
                      len(cluster.qnames), cluster.concordant_frags,
                      ','.join(cluster.qnames)]
            fields = map(str, fields)
            print >>outfh, '\t'.join(fields)
        for r in locus_reads:
            outbamfh.write(r)
    outfh.close()
    outbamfh.close()
    concordant_bamfh.close()
    discordant_bamfh.close()
    #
    # index the newly annotated discordant bam file 
    #
    logging.debug("Indexing newly annotated discordant BAM file")
    pysam.index(output_bam_file)
    #
    # sort the BAM file that has cluster annotations by read name
    #
    logging.debug("Sorting newly annotated discordant BAM file by read name")
    qname_sorted_bam_prefix = os.path.join(tmp_dir, os.path.splitext(output_bam_file)[0] + ".byname")
    qname_sorted_bam_file = qname_sorted_bam_prefix + ".bam"
    pysam.sort("-n", "-m", str(int(1e9)), output_bam_file, qname_sorted_bam_prefix)
    #
    # iterate through named-sorted bam file and aggregate cluster pairs
    #
    logging.debug("Enumerating cluster pairs")
    cluster_pairs = collections.defaultdict(lambda: [])
    bamfh = pysam.Samfile(qname_sorted_bam_file, "rb")
    for pe_reads in parse_pe_reads(bamfh):
        # group into 5' and 3' reads
        reads5p = []
        reads3p = []
        for reads in pe_reads:
            for r in reads:
                orientation = r.opt(ORIENTATION_TAG)
                if orientation == ORIENTATION_5P:
                    reads5p.append(r)
                else:
                    reads3p.append(r)
        # iterate through possible pairs
        for r5p in reads5p:
            for r3p in reads3p:
                cluster_id_5p = r5p.opt(DISCORDANT_CLUSTER_TAG)
                cluster_id_3p = r3p.opt(DISCORDANT_CLUSTER_TAG)
                cluster_pairs[(cluster_id_5p, cluster_id_3p)].append(r5p.qname)
    bamfh.close()
    # write cluster pairs
    outfh = open(cluster_pair_file, "w")
    for idtuple, qnames in cluster_pairs.iteritems():
        id5p, id3p = idtuple
        print >>outfh, '\t'.join(map(str, [id5p, id3p, ','.join(qnames)]))
    outfh.close()
    # remove temporary files
    if os.path.exists(qname_sorted_bam_file):
        os.remove(qname_sorted_bam_file)
    return config.JOB_SUCCESS
    
def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("--tmp-dir", dest="tmp_dir", default=None)
    parser.add_argument("discordant_bam_file") 
    parser.add_argument("concordant_bam_file") 
    parser.add_argument("output_bam_file") 
    parser.add_argument("cluster_file") 
    parser.add_argument("cluster_pair_file") 
    args = parser.parse_args()
    return cluster_discordant_reads(args.discordant_bam_file,
                                    args.concordant_bam_file, 
                                    args.output_bam_file, 
                                    args.cluster_file,
                                    args.cluster_pair_file,
                                    tmp_dir=args.tmp_dir)

if __name__ == '__main__':
    sys.exit(main())