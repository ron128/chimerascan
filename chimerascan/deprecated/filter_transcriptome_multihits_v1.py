'''
Created on Jun 8, 2012

@author: mkiyer
'''
import sys
import logging
import argparse
import os

# local imports
import chimerascan.pysam as pysam
from chimerascan.lib import config
from chimerascan.lib.base import imin2
from chimerascan.lib.sam import parse_pe_reads, copy_read
from chimerascan.lib.feature import TranscriptFeature
from chimerascan.lib.transcriptome_to_genome import \
    build_tid_transcript_genome_map, \
    transcript_to_genome_pos

def annotate_multihits(reads, tid_tx_genome_map):
    hits = set()
    for r in reads:
        if r.is_unmapped:
            continue
        assert r.rname in tid_tx_genome_map
        # use the position that is most 5' relative to genome
        left_tid, left_strand, left_pos = transcript_to_genome_pos(r.rname, r.pos, tid_tx_genome_map)
        right_tid, right_strand, right_pos = transcript_to_genome_pos(r.rname, r.aend-1, tid_tx_genome_map)
        tid = left_tid
        pos = imin2(left_pos, right_pos)
        hits.add((tid, pos))
    for i,r in enumerate(reads):
        # annotate reads with 'HI', and 'IH' tags
        r.tags = r.tags + [("HI",i), ("IH",len(reads)), ("NH",len(hits))]
    return len(hits)

def filter_multihits(transcript_file, input_bam_file, output_bam_file,
                     max_multihits=1):
    logging.debug("Reading transcript features")
    transcripts = list(TranscriptFeature.parse(open(transcript_file)))
    # parse and convert sam -> bam
    inbamfh = pysam.Samfile(input_bam_file, "rb")
    outbamfh = pysam.Samfile(output_bam_file, "wb", template=inbamfh)
    # build a transcript to genome coordinate map   
    tid_tx_genome_map = build_tid_transcript_genome_map(outbamfh, transcripts)
    num_frags = 0
    logging.debug("Annotating and filtering multihits")
    for pe_reads in parse_pe_reads(inbamfh):        
        mate_num_hits = []
        for reads in pe_reads:
            num_hits = annotate_multihits(reads, tid_tx_genome_map)
            mate_num_hits.append(num_hits)
        new_pe_reads = [[],[]]
        if mate_num_hits[0] > max_multihits:
            r = copy_read(pe_reads[0][0])
            r.is_unmapped = True
            r.is_proper_pair = False
            r.is_secondary = False
            r.rname = -1
            r.pos = 0
            if mate_num_hits[1] > max_multihits:
                r.mate_is_unmapped = True
                r.mrnm = -1
                r.mpos = 0
            new_pe_reads[0] = [r]
        else:
            new_pe_reads[0] = pe_reads[0]
        if mate_num_hits[1] > max_multihits:
            r = copy_read(pe_reads[1][0])
            r.is_unmapped = True
            r.is_proper_pair = False
            r.is_secondary = False
            r.rname = -1
            r.pos = 0
            if mate_num_hits[0] > max_multihits:
                r.mate_is_unmapped = True
                r.mrnm = -1
                r.mpos = 0
            new_pe_reads[1] = [r]
        else:
            new_pe_reads[1] = pe_reads[1]
        for reads in pe_reads:
            for r in reads:
                outbamfh.write(r)
        num_frags += 1
    logging.debug("Found %d fragments" % (num_frags))
    inbamfh.close()
    outbamfh.close()
    return config.JOB_SUCCESS

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("--max-multihits", dest="max_multihits", type=int)
    parser.add_argument("transcript_file")
    parser.add_argument("input_bam_file")
    parser.add_argument("output_bam_file") 
    args = parser.parse_args()
    return filter_multihits(args.transcript_file, 
                            args.input_bam_file, 
                            args.output_bam_file,
                            args.max_multihits)

if __name__ == '__main__':
    sys.exit(main())


