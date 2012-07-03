'''
Created on Jul 3, 2012

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
from chimerascan.lib.batch_sort import batch_sort
from chimerascan.lib.sam import parse_pe_reads
from chimerascan.lib.chimera import ORIENTATION_TAG, ORIENTATION_5P, \
    ORIENTATION_3P, DISCORDANT_CLUSTER_TAG

def parse_and_group_cluster_pairs(fh):
    prev_id_5p, prev_id_3p = None,None
    qnames = []
    for line in fh:
        fields = line.strip().split('\t')        
        id5p = fields[0]
        id3p = fields[1]
        if (id5p, id3p) != (prev_id_5p, prev_id_3p):
            if len(qnames) > 0:
                yield prev_id_5p, prev_id_3p, qnames
                qnames = []
            prev_id_5p, prev_id_3p = id5p, id3p
        qnames.append(fields[2])
    if len(qnames) > 0:
        yield id5p, id3p, qnames 

def pair_discordant_clusters(discordant_bam_file, cluster_pair_file, tmp_dir):
    #
    # sort the BAM file that has cluster annotations by read name
    #
    logging.debug("Sorting newly annotated discordant BAM file by read name")
    qname_sorted_bam_prefix = os.path.join(tmp_dir, os.path.splitext(discordant_bam_file)[0] + ".byname")
    qname_sorted_bam_file = qname_sorted_bam_prefix + ".bam"
    pysam.sort("-n", "-m", str(int(1e9)), discordant_bam_file, qname_sorted_bam_prefix)
    #
    # iterate through named-sorted bam file write cluster pairs
    #
    logging.debug("Enumerating cluster pairs")
    tmp_cluster_file = os.path.join(tmp_dir, "tmp_clusters.txt")
    tmp_cluster_fh = open(tmp_cluster_file, 'w')
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
                id5p = r5p.opt(DISCORDANT_CLUSTER_TAG)
                id3p = r3p.opt(DISCORDANT_CLUSTER_TAG)
                print >>tmp_cluster_fh, '\t'.join(map(str, (id5p, id3p, r5p.qname)))
    bamfh.close()
    tmp_cluster_fh.close()
    #
    # sort cluster pairs
    #
    logging.debug("Sorting cluster pairs")
    tmp_sorted_cluster_file = os.path.join(tmp_dir, "tmp_clusters.srt.txt")
    def sortfunc(line):
        fields = line.strip().split('\t')
        return (fields[0], fields[1])
    batch_sort(input=tmp_cluster_file,
               output=tmp_sorted_cluster_file,
               key=sortfunc,
               buffer_size=32000,
               tempdirs=[tmp_dir])
    #
    # write cluster pairs
    #
    logging.debug("Grouping cluster pairs")
    pair_id = 0
    outfh = open(cluster_pair_file, "w")
    for id5p, id3p, qnames in parse_and_group_cluster_pairs(open(tmp_sorted_cluster_file)):
        print >>outfh, '\t'.join(map(str, [pair_id, id5p, id3p, ','.join(qnames)]))
        pair_id += 1
    outfh.close()
    # remove temporary files
    if os.path.exists(qname_sorted_bam_file):
        os.remove(qname_sorted_bam_file)
    if os.path.exists(tmp_cluster_file):
        os.remove(tmp_cluster_file)
    if os.path.exists(tmp_sorted_cluster_file):
        os.remove(tmp_sorted_cluster_file)
    return config.JOB_SUCCESS    

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("--tmp-dir", dest="tmp_dir", default=None)
    parser.add_argument("discordant_bam_file") 
    parser.add_argument("cluster_pair_file") 
    args = parser.parse_args()
    return pair_discordant_clusters(args.discordant_bam_file, 
                                    args.cluster_pair_file, 
                                    args.tmp_dir)

if __name__ == '__main__':
    sys.exit(main())
