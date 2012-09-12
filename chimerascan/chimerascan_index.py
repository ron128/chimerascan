#!/usr/bin/env python
'''
Created on Jan 5, 2011

@author: mkiyer

chimerascan: chimeric transcript discovery using RNA-seq

Copyright (C) 2011-2012 Matthew Iyer

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
import argparse
import logging
import os
import shutil
import subprocess
import sys
import collections
import operator

import pysam

# local imports
from chimerascan.lib.feature import TranscriptFeature
from chimerascan.lib.seq import DNA_reverse_complement
from chimerascan.lib.base import up_to_date, check_executable
from chimerascan.lib import config

BASES_PER_LINE = 50

def split_seq(seq, chars_per_line):
    pos = 0
    newseq = []
    while pos < len(seq):
        if pos + chars_per_line > len(seq):        
            endpos = len(seq)
        else:
            endpos = pos + chars_per_line
        newseq.append(seq[pos:endpos])
        pos = endpos
    return '\n'.join(newseq)

def transcript_features_to_fasta(transcript_feature_file, reference_seq_file):
    ref_fa = pysam.Fastafile(reference_seq_file)
    total = 0
    used = 0
    for g in TranscriptFeature.parse(open(transcript_feature_file)):
        total += 1
        exon_seqs = []
        error_occurred = False
        for start, end in g.exons:
            seq = ref_fa.fetch(g.chrom, start, end)
            if (not seq) or (len(seq) < (end - start)):
                logging.warning("transcript id %d exon %s:%d-%d not found in reference" % 
                                (g.tx_id, g.chrom, start, end))
                error_occurred = True
                break
            exon_seqs.append(seq)
        if error_occurred:
            continue
        used += 1
        # make fasta record
        seq = ''.join(exon_seqs)
        # look for sequences containing only 'N's
        base_counts = collections.Counter(seq)
        valid_bases = sum(base_counts[x] for x in 
                          ("A","T","G","C","a","t","g","c"))
        if valid_bases == 0:
            logging.warning("transcript %d at pos %s:%d-%d lacks valid bases" %
                            (g.tx_id, g.chrom, g.tx_start, g.tx_end))
            continue
        # reverse complement negative stranded sequences
        if g.strand == '-':
            seq = DNA_reverse_complement(seq)
        # break seq onto multiple lines
        seqlines = split_seq(seq, BASES_PER_LINE)
        fa_record = (">%d range=%s:%d-%d strand=%s\n%s" % 
                     (g.tx_id, g.chrom, g.tx_start, g.tx_end, g.strand, 
                      seqlines))
        yield g, fa_record
    logging.info("Used %d/%d gene features" % (used,total))
    ref_fa.close()
    
def find_maximum_feature_overlap(features):
    boundaries = []
    for f in features:
        for start,end in f.exons:
            boundaries.append((start, 1))
            boundaries.append((end, -1))
    boundaries.sort(key=operator.itemgetter(0))
    overlap = 0
    max_overlap = 0
    for b,btype in boundaries:
        overlap += btype
        max_overlap = max(max_overlap, overlap)
    return max_overlap

def create_chimerascan_index(output_dir, 
                             genome_fasta_file, 
                             transcript_feature_file):
    # create output dir if it does not exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        logging.info("Created index directory: %s" % (output_dir))
    # copy reference fasta file to output dir and index it
    dst_genome_fasta_file = os.path.join(output_dir, 
                                         config.GENOME_FASTA_FILE)
    msg = "Adding reference genome"
    if (up_to_date(dst_genome_fasta_file, genome_fasta_file)):
        logging.info("[SKIPPED] %s" % (msg))
    else:
        logging.info(msg)
        shutil.copyfile(genome_fasta_file, dst_genome_fasta_file)
        # index the genome fasta file
        logging.info("Indexing FASTA file")
        fh = pysam.Fastafile(dst_genome_fasta_file)
        fh.close()
    # add gene sequences to index
    dst_transcript_feature_file = os.path.join(output_dir, config.TRANSCRIPT_FEATURE_FILE)
    transcript_fasta_file = os.path.join(output_dir, config.TRANSCRIPTOME_FASTA_FILE)
    multimapping_file = os.path.join(output_dir, config.MAX_MULTIMAPPING_FILE)
    msg = "Building transcriptome sequences and gene features"
    if (up_to_date(dst_transcript_feature_file, transcript_feature_file) and
        up_to_date(transcript_fasta_file, dst_transcript_feature_file) and
        up_to_date(multimapping_file, transcript_feature_file)):
        logging.info("[SKIPPED] %s" % (msg))
    else:
        logging.info(msg)
        # write sequences from gene feature file
        logging.info("Adding transcript sequences")
        fasta_fh = open(transcript_fasta_file, "w")
        tx_fh = open(dst_transcript_feature_file, "w")
        chrom_transcript_dict = collections.defaultdict(lambda: [])
        for t, fa_record in transcript_features_to_fasta(transcript_feature_file, 
                                                         dst_genome_fasta_file):
            print >>tx_fh, str(t)
            print >>fasta_fh, fa_record
            chrom_transcript_dict[t.chrom].append(t)
        tx_fh.close()
        fasta_fh.close()
        # find maximum transcript overlap as this informs alignment 
        # parameters controlling multi-mapping read handling
        max_overlap = 0
        for chrom, transcripts in chrom_transcript_dict.iteritems():            
            overlap = find_maximum_feature_overlap(transcripts)
            max_overlap = max(max_overlap, overlap)
        logging.info("Maximum transcript overlap is %d" % (max_overlap))
        fh = open(multimapping_file, "w")
        print >>fh, max_overlap
        fh.close()
        # index the transcript fasta file
        logging.info("Indexing the Transcriptome FASTA file")
        fh = pysam.Fastafile(transcript_fasta_file)
        fh.close()
    #
    # Build Transcriptome alignment index
    #
    skip = True
    index_files = (os.path.join(output_dir, f) for f in config.TRANSCRIPTOME_BOWTIE2_FILES)
    for f in index_files:
        skip = skip and up_to_date(f, transcript_fasta_file) 
    msg = "Building transcriptome index"
    if skip:
        logging.info("[SKIPPED] %s" % (msg))
    else:
        logging.info(msg)
        bowtie_index_name = os.path.join(output_dir, config.TRANSCRIPTOME_INDEX)
        args = [config.BOWTIE2_BUILD_BIN, transcript_fasta_file, bowtie_index_name]
        if subprocess.call(args) != os.EX_OK:
            logging.error("Failed to create alignment index")
            for f in index_files:
                if os.path.exists(f):
                    os.remove(f)
            return config.JOB_ERROR
    #
    # Build Genome alignment index
    #
    skip = True
    index_files = (os.path.join(output_dir, f) for f in config.GENOME_BOWTIE2_FILES)
    for f in index_files:
        skip = skip and up_to_date(f, dst_genome_fasta_file) 
    msg = "Building genome index"
    if skip:
        logging.info("[SKIPPED] %s" % (msg))
    else:
        logging.info(msg)
        bowtie_index_name = os.path.join(output_dir, config.GENOME_INDEX)
        args = [config.BOWTIE2_BUILD_BIN, dst_genome_fasta_file, bowtie_index_name]
        if subprocess.call(args) != os.EX_OK:
            logging.error("Failed to create alignment index")
            for f in index_files:
                if os.path.exists(f):
                    os.remove(f)
            return config.JOB_ERROR
    logging.info("Chimerascan index created successfully")
    return config.JOB_SUCCESS

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser(description="Build alignment indexes for use with chimerascan")
    parser.add_argument("ref_fasta_file", help="reference genome FASTA file")
    parser.add_argument("transcript_feature_file", help="transcript features")
    parser.add_argument("output_dir", help="directory where indexes will be created")
    args = parser.parse_args()
    # check that input files exist
    if not os.path.isfile(args.ref_fasta_file):
        parser.error("Reference fasta file '%s' not found" % (args.ref_fasta_file))
    if not os.path.isfile(args.transcript_feature_file):
        parser.error("Gene feature file '%s' not found" % (args.transcript_feature_file))
    # check that output dir is not a regular file
    if os.path.exists(args.output_dir) and (not os.path.isdir(args.output_dir)):
        parser.error("Output directory name '%s' exists and is not a valid "
                     "directory" % (args.output_dir))
    # check that bowtie2-build program exists
    if check_executable(config.BOWTIE2_BUILD_BIN):
        logging.debug("Checking for '%s' binary... found" % (config.BOWTIE2_BUILD_BIN))
    else:
        parser.error("%s binary not found or not executable" % (config.BOWTIE2_BUILD_BIN))
    # run main index creation function
    retcode = create_chimerascan_index(args.output_dir, args.ref_fasta_file, 
                                       args.transcript_feature_file)
    return retcode

if __name__ == '__main__':
    sys.exit(main())