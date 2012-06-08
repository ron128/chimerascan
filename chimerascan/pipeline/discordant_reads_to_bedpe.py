'''
Created on Jul 21, 2011

@author: mkiyer
'''
import logging
import os
import sys
import argparse

from chimerascan import pysam
from chimerascan.lib import config
from chimerascan.lib.feature import TranscriptFeature
from chimerascan.lib.chimera import DiscordantTags, DISCORDANT_TAG_NAME, \
    OrientationTags, ORIENTATION_TAG_NAME, DiscordantRead
from chimerascan.lib.transcriptome_to_genome import build_tid_transcript_map
from chimerascan.lib.batch_sort import batch_sort

def parse_pairs(bamfh):
    bam_iter = iter(bamfh)
    try:
        while True:
            r1 = bam_iter.next()
            r2 = bam_iter.next()
            yield r1,r2
    except StopIteration:
        pass
    
def parse_gene_discordant_reads(bamfh):
    """
    return tuples of (5',3') reads that both align to transcripts
    """
    for r1,r2 in parse_pairs(bamfh):
        # TODO:
        # for now we are only going to deal with gene-gene
        # chimeras and leave other chimeras for study at a 
        # later time
        dr1 = r1.opt(DISCORDANT_TAG_NAME)
        dr2 = r2.opt(DISCORDANT_TAG_NAME)
        if (dr1 != DiscordantTags.DISCORDANT_GENE or
            dr2 != DiscordantTags.DISCORDANT_GENE):            
            continue
        # organize key in 5' to 3' order
        or1 = r1.opt(ORIENTATION_TAG_NAME)
        or2 = r2.opt(ORIENTATION_TAG_NAME)
        # TODO: remove assert
        assert or1 != or2
        if or1 == OrientationTags.FIVEPRIME:
            pair = (r1,r2)
        else:
            pair = (r2,r1)
        yield pair

def discordant_reads_to_bedpe(index_dir, input_bam_file, output_file):
    # open BAM alignment file
    bamfh = pysam.Samfile(input_bam_file, "rb")
    # build a lookup table to get genomic intervals from transcripts
    logging.debug("Reading transcript features")
    transcript_file = os.path.join(index_dir, config.TRANSCRIPT_FEATURE_FILE)
    transcripts = list(TranscriptFeature.parse(open(transcript_file)))
    tid_tx_map = build_tid_transcript_map(bamfh, transcripts)
    outfh = open(output_file, "w")    
    logging.debug("Converting BAM to BEDPE format")
    for r5p,r3p in parse_gene_discordant_reads(bamfh):
        # store pertinent read information in lightweight structure called
        # DiscordantRead object. this departs from SAM format into a 
        # custom read format
        dr5p = DiscordantRead.from_read(r5p)
        dr3p = DiscordantRead.from_read(r3p)
        # get gene information
        tx5p = tid_tx_map[r5p.rname]
        tx3p = tid_tx_map[r3p.rname]
        # write bedpe format
        fields = [tx5p.tx_id, r5p.pos, r5p.aend,
                  tx3p.tx_id, r3p.pos, r3p.aend,
                  r5p.qname,  # read name
                  0, # score
                  tx5p.strand, tx3p.strand, # strand 1, strand 2
                  ]
        fields.append('|'.join(map(str, dr5p.to_list())))
        fields.append('|'.join(map(str, dr3p.to_list())))  
        print >>outfh, '\t'.join(map(str, fields)) 
    outfh.close()
    bamfh.close()
    
def sort_bedpe(input_file, output_file, tmp_dir):
    # sort BEDPE file by paired chromosome/position
    def sortfunc(line):
        fields = line.strip().split('\t')
        return tuple([fields[0], fields[3], fields[1], fields[4]])
    tempdirs = [tmp_dir]
    batch_sort(input=input_file,
               output=output_file,
               key=sortfunc,
               buffer_size=32000,
               tempdirs=tempdirs)


def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("index_dir")
    parser.add_argument("bam_file")
    parser.add_argument("bedpe_output_file")
    args = parser.parse_args()
    return discordant_reads_to_bedpe(args.index_dir, 
                                     args.bam_file, 
                                     args.bedpe_output_file) 

if __name__ == '__main__':
    sys.exit(main())