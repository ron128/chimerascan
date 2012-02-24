#!/usr/bin/env python
'''
Created on Jan 5, 2011

@author: mkiyer

chimerascan: chimeric transcript discovery using RNA-seq

Copyright (C) 2011 Matthew Iyer

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
import collections
import logging
import os
import shutil
import subprocess
import sys
from optparse import OptionParser

# local imports
import chimerascan.pysam as pysam
from chimerascan.lib.feature import GeneFeature
from chimerascan.lib.seq import DNA_reverse_complement
from chimerascan.lib.base import up_to_date, check_executable
from chimerascan.bx.intersection import Interval, IntervalTree
from chimerascan.lib.config import JOB_ERROR, JOB_SUCCESS, ALIGN_INDEX, \
    BOWTIE_INDEX_FILE, FRAG_SIZE_INDEX, FRAG_SIZE_INDEX_FILE, \
    GENE_FEATURE_FILE, GENE_REF_PREFIX, RAW_JUNCS_FILE

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

def bed12_to_fasta(gene_feature_file, reference_seq_file):
    ref_fa = pysam.Fastafile(reference_seq_file)
    for g in GeneFeature.parse(open(gene_feature_file)):
        exon_seqs = []
        error_occurred = False
        for start, end in g.exons:
            seq = ref_fa.fetch(g.chrom, start, end)
            if not seq:
                logging.warning("gene %s exon %s:%d-%d not found in reference" % 
                                (g.tx_name, g.chrom, start, end))
                error_occurred = True
                break
            exon_seqs.append(seq)
        if error_occurred:
            continue
        # make fasta record
        seq = ''.join(exon_seqs)
        if g.strand == '-':
            seq = DNA_reverse_complement(seq)
        # break seq onto multiple lines
        seqlines = split_seq(seq, BASES_PER_LINE)    
        yield (">%s range=%s:%d-%d gene=%s strand=%s\n%s" % 
               (GENE_REF_PREFIX + g.tx_name, g.chrom, start, end, g.gene_name, g.strand, seqlines))
    ref_fa.close()

def build_exon_trees(genes):
    trees = collections.defaultdict(lambda: IntervalTree())
    for g in genes:        
        for e in g.exons:
            start, end = e
            trees[g.chrom].insert_interval(Interval(start, end, strand=g.strand))
    return trees

def find_unambiguous_exon_intervals(genes):
    """
    returns (chrom, start, end, strand) tuples for exon
    intervals that are unique and have no overlapping
    transcripts or exons.    
    """
    trees = build_exon_trees(genes)    
    for g in genes:
        for start,end in g.exons:
            hits = [(hit.start, hit.end, hit.strand) 
                    for hit in trees[g.chrom].find(start, end)]
            overlapping_hits = set([(start, end, g.strand)]).union(hits)
            if len(overlapping_hits) == 1:
                yield g.chrom, start, end, g.strand

def create_fragment_size_index(output_dir, gene_feature_file, 
                               reference_seq_file, bowtie_build_bin, 
                               max_fragment_size):
    """
    make an alignment index containing sequences that can be used to
    assess the fragment size distribution.  these sequences must be 
    larger than the 'max_insert_size' in order to be viable for use 
    in characterizing the fragment size distribution.
    """
    # parse genes file
    genes = [g for g in GeneFeature.parse(open(gene_feature_file))]
    # find all exons that are larger than the maximum estimated fragment size
    exons = set([coord for coord in find_unambiguous_exon_intervals(genes)
                 if (coord[2] - coord[1]) >= max_fragment_size])
    logging.info("Found %d exons larger than %d" % (len(exons), max_fragment_size))    
    # extract the nucleotide sequence of the exons
    logging.info("Extracting sequences to use for estimating the fragment "
                 " size distribution")
    ref_fa = pysam.Fastafile(reference_seq_file)    
    frag_size_fa_file = os.path.join(output_dir, "frag_size_seq.fa")
    fh = open(frag_size_fa_file, 'w')
    for chrom, start, end, strand in exons:
        seq = ref_fa.fetch(chrom, start, end)
        if not seq:
            logging.warning("exon %s:%d-%d not found in reference" % (chrom, start, end))
            continue
        # make fasta record
        if strand == '-':
            seq = DNA_reverse_complement(seq)
            # break seq onto multiple lines
            seqlines = split_seq(seq, BASES_PER_LINE)    
            record = (">%s:%d-%d strand=%s\n%s" % 
                      (chrom, start, end, strand, seqlines))
            print >>fh, record
    fh.close()
    ref_fa.close()
    # build bowtie alignment index from the fragment size exons
    logging.info("Building bowtie index")
    frag_size_index = os.path.join(output_dir, FRAG_SIZE_INDEX)
    args = [bowtie_build_bin, frag_size_fa_file, frag_size_index]
    return subprocess.call(args)

def create_tophat_juncs_file(output_dir, gene_feature_file):
    """
    adapted from the 'bed_to_juncs' script distributed with the
    TopHat package. http://tophat.cbcb.umd.edu
    """
    line_num = 0
    for line in open(gene_feature_file):
        line = line.strip()
        if line.startswith("#"):
            continue
        fields = line.split()
        if len(fields) < 10:
            logging.warning("Malformed line %d, missing columns" % (line_num))
            continue
        line_num += 1
        chrom = fields[1]
        strand = fields[2]
        tx_start = int(fields[3])
        #tx_end = int(fields[4])
        exon_starts = map(int, fields[8].split(",")[:-1])
        exon_ends = map(int, fields[9].split(",")[:-1])
        for i in xrange(1,len(exon_starts)):
            junc_start = tx_start + exon_ends[i-1] - 1
            junc_end = tx_start + exon_starts[i]
            yield "%s\t%d\t%d\t%s" % (chrom, junc_start, junc_end, strand)

def create_chimerascan_index(output_dir, 
                             genome_fasta_file, 
                             gene_feature_file,
                             bowtie_build_bin):
#                             min_fragment_size,
#                             max_fragment_size):
    # create output dir if it does not exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        logging.info("Created index directory: %s" % (output_dir))
    # copy reference fasta file to output dir
    index_fasta_file = os.path.join(output_dir, ALIGN_INDEX + ".fa")
    if (up_to_date(index_fasta_file, genome_fasta_file) and
        up_to_date(index_fasta_file, gene_feature_file)):
        logging.info("[SKIPPED] Adding reference genome to index")
    else:
        logging.info("Adding reference genome to index")
        shutil.copyfile(genome_fasta_file, index_fasta_file)
        # index the genome fasta file
        logging.info("Indexing FASTA file")
        fh = pysam.Fastafile(index_fasta_file)
        fh.close()
        # append sequences from gene feature file
        logging.info("Adding transcript sequences to index...")
        fh = open(index_fasta_file, "a")
        for fa_record in bed12_to_fasta(gene_feature_file, 
                                        index_fasta_file):
            print >>fh, fa_record
        fh.close()
        # remove old fasta index
        os.remove(index_fasta_file + ".fai")
        # re-index the combined fasta file
        logging.info("Re-indexing FASTA file...")
        fh = pysam.Fastafile(index_fasta_file)
        fh.close()
    # build bowtie index on the reference sequence file
    bowtie_index_file = os.path.join(output_dir, BOWTIE_INDEX_FILE)
    msg = "Building bowtie index"
    if up_to_date(bowtie_index_file, index_fasta_file):
        logging.info("[SKIPPED] %s" % (msg))
    else:
        logging.info(msg)
        bowtie_index_name = os.path.join(output_dir, ALIGN_INDEX)
        args = [bowtie_build_bin, index_fasta_file, bowtie_index_name]
        if subprocess.call(args) != os.EX_OK:
            logging.error("bowtie-build failed to create alignment index")
            if os.path.exists(bowtie_index_file):
                os.remove(bowtie_index_file)
            return JOB_ERROR
    # copy gene bed file to index directory
    dst_gene_feature_file = os.path.join(output_dir, GENE_FEATURE_FILE)
    if up_to_date(dst_gene_feature_file, gene_feature_file):
        logging.info("[SKIPPED] Adding transcript features to index...")
    else:
        logging.info("Adding transcript features to index...")
        shutil.copyfile(gene_feature_file, dst_gene_feature_file)
    # create tophat junctions file from gene features
#    juncs_file = os.path.join(output_dir, TOPHAT_JUNCS_FILE)
#    if up_to_date(juncs_file, dst_gene_feature_file):
#        logging.info("[SKIPPED] Creating splice junction file...")
#    else:
#        logging.info("Creating splice junction file...")
#        fh = open(juncs_file, "w")
#        for junc_line in create_tophat_juncs_file(output_dir, gene_feature_file):
#            print >>fh, junc_line
#        fh.close()
    # build special index used to discover the fragment size
#    frag_size_index_file = os.path.join(output_dir, FRAG_SIZE_INDEX_FILE)
#    if up_to_date(frag_size_index_file, index_fasta_file):
#        logging.info("[SKIPPED] Building fragment size distribution index")
#    else:
#        logging.info("Building fragment size distribution index")
#        retcode = create_fragment_size_index(output_dir, gene_feature_file, 
#                                             genome_fasta_file, 
#                                             bowtie_build_bin, 
#                                             max_fragment_size)
#        if retcode != os.EX_OK:
#            logging.error("bowtie-build failed to create fragment size "
#                          "distribution index")
#            if os.path.exists(frag_size_index_file):
#                os.remove(frag_size_index_file)
#            return JOB_ERROR 
    logging.info("chimerascan index created successfully")
    return JOB_SUCCESS


def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = OptionParser("usage: %prog [options] <reference_genome.fa> "
                          "<gene_models.txt> <index_output_dir>")
    #parser.add_option('-i', '--min-fragment-size', dest="min_fragment_size", default=0)
    #parser.add_option('-I', '--max-fragment-size', dest="max_fragment_size", default=700)    
    parser.add_option("--bowtie-dir", dest="bowtie_dir", default="",
                      help="Path to the 'bowtie' software (by default, "
                      "expects the 'bowtie' and 'bowtie-build' "
                      "binaries to be in current PATH)")
    options, args = parser.parse_args()
    # check command line arguments
    if len(args) < 3:
        parser.error("Incorrect number of command line arguments")
    ref_fasta_file = args[0]
    gene_feature_file = args[1]
    output_dir = args[2]
    # check that input files exist
    if not os.path.isfile(ref_fasta_file):
        parser.error("Reference fasta file '%s' not found" % (ref_fasta_file))
    if not os.path.isfile(gene_feature_file):
        parser.error("Gene feature file '%s' not found" % (gene_feature_file))
    # check that output dir is not a regular file
    if os.path.exists(output_dir) and (not os.path.isdir(output_dir)):
        parser.error("Output directory name '%s' exists and is not a valid "
                     "directory" % (output_dir))
    # check that bowtie-build program exists
    bowtie_build_bin = os.path.join(options.bowtie_dir, "bowtie-build")
    if check_executable(bowtie_build_bin):
        logging.debug("Checking for 'bowtie-build' binary... found")
    else:
        parser.error("bowtie-build binary not found or not executable")
    # run main index creation function
    retcode = create_chimerascan_index(output_dir, ref_fasta_file, 
                                       gene_feature_file,
                                       bowtie_build_bin)
#                                       min_fragment_size=options.min_fragment_size,
#                                       max_fragment_size=options.max_fragment_size)
    sys.exit(retcode)


if __name__ == '__main__':
    main()