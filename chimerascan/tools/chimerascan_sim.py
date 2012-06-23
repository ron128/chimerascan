'''
Created on Feb 24, 2012

@author: mkiyer
'''
import os
import sys
import argparse
import logging
import operator
import collections
import random

import pysam

from chimerascan.lib import config
from chimerascan.lib.seq import DNA_reverse_complement
from chimerascan.lib.feature import TranscriptFeature

DEFAULT_FRAG_SIZE_MEAN = 200
DEFAULT_FRAG_SIZE_SD = 20
DEFAULT_READ_LENGTH = 50
DEFAULT_NUM_FRAGS = 1000000
DEFAULT_STRANDED = False

def parse_transcript_exprs_file(fh,
                                transcript_id_field="tracking_id",
                                fpkm_field="FPKM"):
    header_fields = fh.next().strip().split('\t')
    transcript_id_ind = header_fields.index(transcript_id_field)
    fpkm_ind = header_fields.index(fpkm_field)
    for line in fh:
        fields = line.strip().split('\t')
        transcript_id = fields[transcript_id_ind]
        fpkm = float(fields[fpkm_ind])
        yield transcript_id, fpkm

def get_transcript_sequence(fastafh, transcript):
    exon_seqs = []
    for start,end in transcript.exons:
        seq = fastafh.fetch(transcript.chrom, start, end)
        if (not seq) or (len(seq) < (end - start)):
            logging.warning("exon %s:%d-%d not found in reference" % 
                            (transcript.chrom, start, end))
            return None
        exon_seqs.append(seq)
    seq = ''.join(exon_seqs)
    if transcript.strand == "-":
        seq = DNA_reverse_complement(seq)
    return seq

def randomize_strand(left_pos, left_seq, right_pos, right_seq):
    if random.choice((False,True)):
        return left_pos, left_seq, right_pos, right_seq 
    else:
        return right_pos, right_seq, left_pos, left_seq

def generate_random_frags(seq, frags, frag_size_mean, frag_size_sd, 
                          read_length, num_frags, stranded):
    for i in xrange(frags):
        frag_size = int(round(random.normalvariate(frag_size_mean, frag_size_sd)))
        if frag_size > len(seq):
            continue
        left_pos = random.randint(0, len(seq) - frag_size)
        left_seq = seq[left_pos:left_pos+read_length]
        right_pos = left_pos + frag_size - read_length
        right_seq = seq[right_pos:right_pos+read_length]
        right_seq = DNA_reverse_complement(right_seq)
        if not stranded:
            yield randomize_strand(left_pos, left_seq, right_pos, right_seq)
        else:
            yield left_pos, left_seq, right_pos, right_seq

def generate_transcript_reads(fastafh,
                              transcript_dict, 
                              transcript_exprs_file,
                              frag_size_mean,
                              frag_size_sd,
                              read_length,
                              num_frags,
                              stranded):
    for tx_id, fpkm in \
        parse_transcript_exprs_file(open(transcript_exprs_file)):
        # make full length transcript sequence from exon features 
        seq = get_transcript_sequence(fastafh, transcript_dict[tx_id])
        if seq is None:
            logging.warning("could not extract sequence for transcript %s" % (tx_id))
            continue
        # figure out how many fragments to generate
        frags = int(round(fpkm * (len(seq) / 1000.0) * (num_frags / 1.0e6)))
        for frag_tuple in generate_random_frags(seq, frags, frag_size_mean, 
                                                frag_size_sd, read_length, 
                                                num_frags, stranded):            
            yield (tx_id,) + frag_tuple

def parse_chimera_file(fh):
    for line in fh:
        if line.startswith("#"):
            continue
        fields = line.strip().split('\t')        
        t5p = fields[0]
        start5p = int(fields[1])
        end5p = int(fields[2])
        t3p = fields[3]
        start3p = int(fields[4])
        end3p = int(fields[5])
        fpkm = float(fields[6])
        yield (t5p, start5p, end5p, t3p, start3p, end3p, fpkm)

def generate_fusion_reads(fastafh, transcript_dict, chimera_file,
                          frag_size_mean,
                          frag_size_sd,
                          read_length,
                          num_frags,
                          stranded):
    for t1,s1,e1,t2,s2,e2,fpkm in parse_chimera_file(open(chimera_file)):
        name = "%s:%d-%d|%s:%d-%d" % (t1,s1,e1,t2,s2,e2)
        t1seq = get_transcript_sequence(fastafh, transcript_dict[t1])
        t2seq = get_transcript_sequence(fastafh, transcript_dict[t2])
        if (t1seq is None) or (t2seq is None):
            continue
        t1seq = t1seq[s1:e1]
        t2seq = t2seq[s2:e2]
        seq = t1seq + t2seq
        junc_pos = len(t1seq)
        # figure out how many fragments to generate
        frags = int(round(fpkm * (len(seq) / 1000.0) * (num_frags / 1.0e6)))
        for pos1,seq1,pos2,seq2 in generate_random_frags(seq, frags, frag_size_mean, 
                                                         frag_size_sd, read_length, 
                                                         num_frags, stranded):
            left,right = (pos1,pos2) if (pos1 <= pos2) else (pos2,pos1)
            brk1 = (pos1 < junc_pos) and ((pos1+read_length) > junc_pos)
            brk2 = (pos2 < junc_pos) and ((pos2+read_length) > junc_pos)
            span = (left < junc_pos) and ((right+read_length) > junc_pos)
            yield (name, span, brk1, brk2, pos1, seq1, pos2, seq2)

def to_fastq(qname, readnum, seq, qual):
    return "@%s/%d\n%s\n+\n%s" % (qname, readnum, seq, qual)
    
def generate_simulated_reads(fastafh,
                             transcript_dict,
                             transcript_exprs_file, 
                             chimera_file, 
                             output_prefix,
                             frag_size_mean=DEFAULT_FRAG_SIZE_MEAN, 
                             frag_size_sd=DEFAULT_FRAG_SIZE_SD,
                             num_frags=DEFAULT_NUM_FRAGS,
                             read_length=DEFAULT_READ_LENGTH,
                             stranded=DEFAULT_STRANDED):
    # generate baseline transcript reads
    logging.info("Generating background reads")
    fragnum = 1
    fh1 = open(output_prefix + "_1.fq", "w")
    fh2 = open(output_prefix + "_2.fq", "w")
    for name, r1pos, r1seq, r2pos, r2seq in \
        generate_transcript_reads(fastafh, transcript_dict,
                                  transcript_exprs_file, frag_size_mean, 
                                  frag_size_sd, read_length, num_frags, 
                                  stranded):
        qual = "I" * read_length
        r1name = "R%d:%s:%d-%d" % (fragnum, name, r1pos, r1pos+read_length)
        r2name = "R%d:%s:%d-%d" % (fragnum, name, r2pos, r2pos+read_length)
        print >>fh1, to_fastq(r1name, 1, r1seq, qual)
        print >>fh2, to_fastq(r2name, 2, r2seq, qual)
        fragnum += 1
    # generate fusion reads
    logging.info("Generating fusion reads")
    for name, span, brk1, brk2, pos1, seq1, pos2, seq2 in \
        generate_fusion_reads(fastafh, transcript_dict, chimera_file,
                              frag_size_mean, frag_size_sd, read_length, 
                              num_frags, stranded):
        qual = "I" * read_length
        r1name = ("F%d:%s:SPAN%d:R1B%d:R2B%d:%d-%d" % 
                  (fragnum, name, int(span), int(brk1), int(brk2), 
                   pos1, pos1+read_length))
        r2name = ("F%d:%s:SPAN%d:R1B%d:R2B%d:%d-%d" % 
                  (fragnum, name, int(span), int(brk1), int(brk2), 
                   pos2, pos2+read_length))
        print >>fh1, to_fastq(r1name, 1, seq1, qual)
        print >>fh2, to_fastq(r2name, 2, seq2, qual)
        fragnum += 1
    fastafh.close()
    fh1.close()
    fh2.close()

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("--frag-size-mean", dest="frag_size_mean", 
                        type=float, default=DEFAULT_FRAG_SIZE_MEAN)
    parser.add_argument("--frag-size-sd", dest="frag_size_sd", 
                        type=float, default=DEFAULT_FRAG_SIZE_SD)
    parser.add_argument("--rlen", dest="read_length", 
                        type=int, default=DEFAULT_READ_LENGTH)
    parser.add_argument("--stranded", dest="stranded", action="store_true",
                        default=False)
    parser.add_argument("-n", dest="num_frags", type=int, 
                        default=DEFAULT_NUM_FRAGS)
    parser.add_argument("index_dir")
    parser.add_argument("transcript_exprs_file")
    parser.add_argument("chimera_file")
    parser.add_argument("output_prefix")
    args = parser.parse_args()
    # setup parameters
    transcript_file = os.path.join(args.index_dir, config.TRANSCRIPT_FEATURE_FILE)
    genome_fasta_file = os.path.join(args.index_dir, config.GENOME_FASTA_FILE)
    fastafh = pysam.Fastafile(genome_fasta_file)
    logging.info("Reading transcripts")
    transcript_dict = {}
    for t in TranscriptFeature.parse(open(transcript_file)):
        transcript_dict[str(t.tx_id)] = t
    logging.info("Generating simulated reads")
    generate_simulated_reads(fastafh,
                             transcript_dict,
                             args.transcript_exprs_file, 
                             args.chimera_file, 
                             args.output_prefix,
                             frag_size_mean=args.frag_size_mean,
                             frag_size_sd=args.frag_size_sd,
                             num_frags=args.num_frags,
                             read_length=args.read_length,
                             stranded=args.stranded)
    fastafh.close()

if __name__ == '__main__':
    sys.exit(main())
