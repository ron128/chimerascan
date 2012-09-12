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

#def _fetch_cluster_boundary_reads(bamfh, qnames, cluster):
#    boundary_pos = _get_cluster_boundary(cluster)
#    reads = []
#    # fetch discordant pair reads in cluster
#    for r in bamfh.fetch(cluster.rname, cluster.start, cluster.end):
#        # read must be part of this cluster-pair
#        if r.qname not in qnames:
#            continue
#        # untrimmed read should overlap cluster boundary
#        padstart, padend = get_clipped_interval(r)
#        if padstart < boundary_pos < padend:
#            reads.append(r)
#    return reads

def _fetch_cluster_reads(bamfh, qnames, cluster):
    reads = []
    # fetch discordant pair reads in cluster
    for r in bamfh.fetch(cluster.rname, cluster.start, cluster.end):
        # read must be part of this cluster-pair
        if r.qname not in qnames:
            continue
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
    reads.extend(_fetch_cluster_reads(discordant_bamfh, qnames, cluster5p))
    reads.extend(_fetch_cluster_reads(discordant_bamfh, qnames, cluster3p))
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


def realign_across_breakpoints(index_dir, 
                               discordant_bam_file,
                               unpaired_bam_file,
                               cluster_shelve_file, 
                               cluster_pair_file, 
                               breakpoint_bam_file,
                               log_dir,
                               tmp_dir,
                               num_processors,
                               local_anchor_length,
                               local_multihits):
    # load cluster database file
    cluster_shelve = shelve.open(cluster_shelve_file, 'r')
    # open discordant reads file
    discordant_bamfh = pysam.Samfile(discordant_bam_file, "rb")
    unpaired_bamfh = pysam.Samfile(unpaired_bam_file, "rb")
    # create tmp dir if it does not exist
    fastq_file = os.path.join(tmp_dir, config.BREAKPOINT_FASTQ_FILE)
    fastq_fh = open(fastq_file, 'w')
    # iterate through cluster pairs and get breakpoint reads
    logging.debug("Extracting breakpoint spanning sequences")
    num_seqs = 0
    for cluster_pair in parse_discordant_cluster_pair_file(open(cluster_pair_file)):
        for fastq_line in _get_cluster_breakpoint_fastq(cluster_pair, 
                                                        cluster_shelve, 
                                                        discordant_bamfh, 
                                                        unpaired_bamfh):
            print >>fastq_fh, fastq_line
            num_seqs += 1
    fastq_fh.close()        
    discordant_bamfh.close()
    unpaired_bamfh.close()
    logging.debug("\tFound %d putative breakpoint spanning sequences" % (num_seqs))
    # use bowtie2 local alignment to find spanning reads 
    transcriptome_index = os.path.join(index_dir, config.TRANSCRIPTOME_INDEX)
    genome_index = os.path.join(index_dir, config.GENOME_INDEX)
    transcript_file = os.path.join(index_dir, config.TRANSCRIPT_FEATURE_FILE)
    log_file = os.path.join(log_dir, config.BREAKPOINT_LOG_FILE)
    logging.debug("Realigning breakpoint spanning sequences")
    bowtie2_align_local(transcriptome_index,
                        genome_index,
                        transcript_file,                                   
                        fastq_file,
                        breakpoint_bam_file,
                        log_file,
                        local_anchor_length=local_anchor_length,
                        local_multihits=local_multihits,
                        num_processors=num_processors)
    cluster_shelve.close()
    return config.JOB_SUCCESS


def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("--tmp-dir", dest="tmp_dir", default="/tmp")
    parser.add_argument("--log-dir", dest="log_dir", default="/tmp")
    parser.add_argument("-p", type=int, dest="num_processors", default=config.BASE_PROCESSORS)
    parser.add_argument("--local-anchor-length", type=int, 
                        dest="local_anchor_length", 
                        default=config.DEFAULT_LOCAL_ANCHOR_LENGTH)
    parser.add_argument("--local-multihits", type=int, 
                        dest="local_multihits", 
                        default=config.DEFAULT_LOCAL_MULTIHITS)
    parser.add_argument("index_dir")
    parser.add_argument("discordant_bam_file")    
    parser.add_argument("unpaired_bam_file")    
    parser.add_argument("cluster_shelve_file")
    parser.add_argument("cluster_pair_file")
    parser.add_argument("breakpoint_bam_file")
    args = parser.parse_args()    
    # run main function
    retcode = realign_across_breakpoints(args.index_dir, 
                                         args.discordant_bam_file,
                                         args.unpaired_bam_file,
                                         args.cluster_shelve_file, 
                                         args.cluster_pair_file, 
                                         args.breakpoint_bam_file,
                                         log_dir=args.log_dir,
                                         tmp_dir=args.tmp_dir,
                                         num_processors=args.num_processors,
                                         local_anchor_length=args.local_anchor_length,
                                         local_multihits=args.local_multihits)
    return retcode

if __name__ == "__main__":
    sys.exit(main())
