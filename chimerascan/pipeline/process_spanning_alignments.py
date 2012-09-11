'''
Created on Sep 10, 2012

@author: mkiyer
'''

'''
Created on Jul 8, 2012

@author: mkiyer
'''
import argparse
import logging
import os
import sys
import shelve

import pysam

from chimerascan.lib import config
from chimerascan.lib.chimera import Chimera, \
    parse_discordant_cluster_pair_file, ORIENTATION_5P, ORIENTATION_3P
from chimerascan.lib.feature import TranscriptFeature
from chimerascan.lib.sam import get_clipped_interval, parse_reads_by_qname
from chimerascan.lib.seq import DNA_reverse_complement
from chimerascan.pipeline.align_bowtie2 import bowtie2_align_local

import chimerascan.pipeline
_pipeline_dir = chimerascan.pipeline.__path__[0]

def _get_cluster_boundary(cluster):
    if (((cluster.strand == "+") and (cluster.orientation == ORIENTATION_5P)) or
        ((cluster.strand == "-") and (cluster.orientation == ORIENTATION_3P))):
        return cluster.end
    else:
        return cluster.start

def _fetch_cluster_boundary_reads(bamfh, qnames, cluster):
    boundary_pos = _get_cluster_boundary(cluster)
    reads = []
    # fetch discordant pair reads in cluster
    for r in bamfh.fetch(cluster.rname, cluster.start, cluster.end):
        # read must be part of this cluster-pair
        if r.qname not in qnames:
            continue
        # untrimmed read should overlap cluster boundary
        padstart, padend = get_clipped_interval(r)
        if padstart < boundary_pos < padend:
            reads.append(r)
    return reads

def _fetch_unpaired_mates(bamfh, cluster):
    # fetch unpaired reads in 5' cluster
    reads = []
    qnames = set(cluster.unpaired_qnames)
    for r in bamfh.fetch(cluster.rname, cluster.start, cluster.end):
        # read must be part of this cluster-pair
        if r.qname not in qnames:
            continue
        reads.append(r)
    return reads

def _get_fastq(qname, rnum, seq, qual):
    return "@%s/%d\n%s\n+\n%s" % (qname, rnum, seq, qual)

def _get_cluster_breakpoint_fastq(cluster_pair, cluster_shelve, 
                                  discordant_bamfh, 
                                  unpaired_bamfh):
    # lookup 5' and 3' clusters
    cluster5p = cluster_shelve[str(cluster_pair.id5p)]
    cluster3p = cluster_shelve[str(cluster_pair.id3p)]
    # find paired reads overlapping edges
    qnames = set(cluster_pair.qnames)
    reads = []
    reads.extend(_fetch_cluster_boundary_reads(discordant_bamfh, qnames, cluster5p))
    reads.extend(_fetch_cluster_boundary_reads(discordant_bamfh, qnames, cluster3p))
    # yield fastq strings
    for r in reads:
        qname = "%d:%s" % (cluster_pair.pair_id, r.qname)
        rnum = int(r.is_read2) + 1
        if r.is_reverse:
            seq = DNA_reverse_complement(r.seq)
            qual = r.qual[::-1]
        else:
            seq = r.seq
            qual = r.qual
        yield _get_fastq(qname, rnum, seq, qual)
    # find reads within cluster with unmapped mates
    reads = []
    reads.extend(_fetch_unpaired_mates(unpaired_bamfh, cluster5p))
    reads.extend(_fetch_unpaired_mates(unpaired_bamfh, cluster3p))
    for r in reads:
        qname = "%d:%s" % (cluster_pair.pair_id, r.qname)
        rnum = 1 if r.is_read2 else 2
        seq = r.opt('R2')
        qual = r.opt('Q2')
        yield _get_fastq(qname, rnum, seq, qual)

def _parse_bam_by_cluster_pair(bamfh):
    reads = []
    current_pair_id = None
    for r in bamfh:
        pair_id, qname = r.qname.split(':')
        pair_id = int(pair_id)
        r.qname = qname
        if (pair_id != current_pair_id) and (len(reads) > 0):
            yield current_pair_id, reads
            reads = []
        current_pair_id = pair_id
        reads.append(r)
    if len(reads) > 0:
        yield current_pair_id, reads
        reads = []

def _test_read_in_cluster(r, rname, cluster):
    if rname != cluster.rname:
        return False
    if (r.pos < cluster.end) and (r.aend > cluster.start):
        return True
    return False
    
def nominate_spanning_reads(cluster_pair, cluster_shelve, bamfh, cluster_reads):
    # lookup 5' and 3' clusters
    cluster5p = cluster_shelve[str(cluster_pair.id5p)]
    cluster3p = cluster_shelve[str(cluster_pair.id3p)]
    # iterate through cluster pair reads
    spanning_reads = []
    for reads in parse_reads_by_qname(cluster_reads):
        if len(reads) < 2:
            continue
        hits5p = []
        hits3p = []
        for r in reads:
            if r.is_unmapped:
                continue
            rname = bamfh.getrname(r.tid)
            is5p = _test_read_in_cluster(r, rname, cluster5p)
            is3p = _test_read_in_cluster(r, rname, cluster3p)
            if is5p and is3p:
                logging.warning("Read %s has local alignments to both 5' and 3' clusters" % (r.qname))
            elif is5p:
                hits5p.append(r)
            elif is3p:
                hits3p.append(r)
            if (len(hits5p) > 0) and (len(hits3p) > 0):
                # pull out best scoring pair of hits
                spanning_reads.append((hits5p[0],hits3p[0]))
                break
    return spanning_reads

def process_spanning_alignments(cluster_shelve_file, 
                                cluster_pair_file,
                                bam_file, 
                                output_sam_file,
                                output_cluster_pair_file):
    # load cluster database file
    cluster_shelve = shelve.open(cluster_shelve_file, 'r')
    # parse breakpoint alignments and output spanning reads
    bamfh = pysam.Samfile(bam_file, "rb")
    outsamfh = pysam.Samfile(output_sam_file, "wh", template=bamfh)
    outfh = open(output_cluster_pair_file, "w")
    cluster_pair_iter = parse_discordant_cluster_pair_file(open(cluster_pair_file))
    # get cluster reads from BAM file
    num_spanning_reads = 0
    for pair_id, cluster_reads in _parse_bam_by_cluster_pair(bamfh):
        # synch with cluster pair file
        cluster_pair = cluster_pair_iter.next()
        while pair_id != cluster_pair.pair_id:
            # no spanning reads here
            print >>outfh, '\t'.join(map(str, [cluster_pair.pair_id, 
                                               cluster_pair.id5p, 
                                               cluster_pair.id3p, 
                                               ','.join(cluster_pair.qnames),
                                               '']))            
            cluster_pair = cluster_pair_iter.next()
        # get spanning read alignments
        spanning_reads = nominate_spanning_reads(cluster_pair, 
                                                 cluster_shelve, 
                                                 bamfh,
                                                 cluster_reads)
        spanning_qnames = sorted(set(r5p.qname for r5p,r3p in spanning_reads))
        # write new cluster pair file
        print >>outfh, '\t'.join(map(str, [cluster_pair.pair_id, 
                                           cluster_pair.id5p, 
                                           cluster_pair.id3p, 
                                           ','.join(cluster_pair.qnames),
                                           ','.join(spanning_qnames)]))
        # write spanning reads to SAM file
        for r5p,r3p in spanning_reads:
            outsamfh.write(r5p)
            outsamfh.write(r3p)
        num_spanning_reads += len(spanning_reads)
    # finish outputting remaining clusters
    for cluster_pair in cluster_pair_iter:
        print >>outfh, '\t'.join(map(str, [cluster_pair.pair_id, 
                                           cluster_pair.id5p, 
                                           cluster_pair.id3p, 
                                           ','.join(cluster_pair.qnames), 
                                           '']))
    logging.debug("\tFound %d spanning read alignments" % (num_spanning_reads))
    outsamfh.close()
    outfh.close()
    bamfh.close()
    cluster_shelve.close()
    return config.JOB_SUCCESS


def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("cluster_shelve_file")
    parser.add_argument("cluster_pair_file")
    parser.add_argument("bam_file")
    parser.add_argument("output_sam_file")
    parser.add_argument("output_cluster_pair_file")
    args = parser.parse_args()    
    # run main function
    retcode = process_spanning_alignments(args.cluster_shelve_file, 
                                          args.cluster_pair_file,
                                          args.bam_file, 
                                          args.output_sam_file,
                                          args.output_cluster_pair_file)
    return retcode

if __name__ == "__main__":
    sys.exit(main())
