'''
Created on Sep 10, 2012

@author: mkiyer
'''
import argparse
import logging
import sys
import shelve
import operator

import pysam

from chimerascan.lib import config
from chimerascan.lib.chimera import parse_discordant_cluster_pair_file, \
    ORIENTATION_5P, ORIENTATION_3P
from chimerascan.lib.sam import get_clipped_interval, \
    parse_reads_by_qname, CIGAR
from chimerascan.lib.seq import DNA_reverse_complement
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

def _get_best_aligned_seq_interval(read):
    start = 0
    end = 0
    best_interval = None
    best_length = 0
    # account for reverse complemented alignments    
    cigar = list(read.cigar)
    if read.is_reverse:
        cigar.reverse()
    for op,bp in read.cigar:
        if ((op == CIGAR.M) or (op == CIGAR.I) or 
            (op == CIGAR.E) or (op == CIGAR.X)):
            end += bp
        elif op == CIGAR.S:
            length = (end - start)
            if (best_interval is None) or (length > best_length):
                best_interval = (start,end)
                best_length = length
            start = end + bp
            end = start
    length = (end - start)
    if (best_interval is None) or (length > best_length):
        best_interval = (start, end)
        best_length = length
    return best_interval

def _test_interval_overlap(start5p, end5p, start3p, end3p,
                           min_interval_length):
    # test if intervals overlap
    if (start5p < end3p) and (end5p > start3p):
        left, right = sorted((start5p, end5p, start3p, end3p))[1:3]
        overlap = right - left
        # check smallest of two intervals
        trim5p = (end5p - start5p) - overlap
        trim3p = (end3p - start3p) - overlap
        if ((trim5p < min_interval_length) or 
            (trim3p < min_interval_length)):
            return True
    return False
        
def _find_compatible_split_reads(hits5p, hits3p, local_anchor_length):
    pairs = []
    for r5p in hits5p:
        # get 5' aligned interval
        interval5p = _get_best_aligned_seq_interval(r5p)
        if interval5p is None:
            continue
        start5p, end5p = interval5p
        score5p = r5p.opt('AS')
        for r3p in hits3p:
            # get 3' aligned interval
            interval3p = _get_best_aligned_seq_interval(r3p)
            if interval3p is None:
                continue
            start3p, end3p = interval3p
            score3p = r3p.opt('AS')
            if _test_interval_overlap(start5p, end5p, start3p, end3p,
                                      local_anchor_length):
                continue
            # TODO: modify the alignment records!
            pairs.append((score5p + score3p, r5p, r3p))
    pairs.sort(key=operator.itemgetter(0), reverse=True)
    return [(r5p,r3p) for (score,r5p,r3p) in pairs]

def nominate_spanning_reads(cluster_pair, cluster_shelve, bamfh, 
                            cluster_reads, local_anchor_length):
    # lookup 5' and 3' clusters
    cluster5p = cluster_shelve[str(cluster_pair.id5p)]
    cluster3p = cluster_shelve[str(cluster_pair.id3p)]
    # iterate through cluster pair reads
    spanning_reads = []
    for reads in parse_reads_by_qname(cluster_reads):
        if len(reads) < 2:
            continue
        # group reads by cluster
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
        # find compatible pairs of split reads
        pairs = _find_compatible_split_reads(hits5p, hits3p, local_anchor_length)
        if len(pairs) > 0:
            spanning_reads.append(pairs[0])
    return spanning_reads

def process_spanning_alignments(cluster_shelve_file, 
                                cluster_pair_file,
                                bam_file, 
                                output_sam_file,
                                output_cluster_pair_file,
                                local_anchor_length):
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
                                                 cluster_reads,
                                                 local_anchor_length)
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
    parser.add_argument("--local-anchor-length", type=int, 
                        dest="local_anchor_length", 
                        default=config.DEFAULT_LOCAL_ANCHOR_LENGTH)
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
                                          args.output_cluster_pair_file,
                                          args.local_anchor_length)
    return retcode

if __name__ == "__main__":
    sys.exit(main())
