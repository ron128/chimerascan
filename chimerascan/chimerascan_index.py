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
from chimerascan.lib.config import JOB_ERROR, JOB_SUCCESS, ALIGN_INDEX, \
    BOWTIE_INDEX_FILE, GENE_FEATURE_FILE, GENE_REF_PREFIX

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

def genepred_to_fasta(gene_feature_file, reference_seq_file):
    ref_fa = pysam.Fastafile(reference_seq_file)
    total = 0
    used = 0
    for g in GeneFeature.parse(open(gene_feature_file)):
        total += 1
        exon_seqs = []
        error_occurred = False
        for start, end in g.exons:
            seq = ref_fa.fetch(g.chrom, start, end)
            if (not seq) or (len(seq) < (end - start)):
                logging.warning("gene %s exon %s:%d-%d not found in reference" % 
                                (g.tx_name, g.chrom, start, end))
                error_occurred = True
                break
            exon_seqs.append(seq)
        if error_occurred:
            continue
        used += 1
        # make fasta record
        seq = ''.join(exon_seqs)
        if g.strand == '-':
            seq = DNA_reverse_complement(seq)
        # break seq onto multiple lines
        seqlines = split_seq(seq, BASES_PER_LINE)    
        fa_record = (">%s range=%s:%d-%d gene=%s strand=%s\n%s" % 
                     (GENE_REF_PREFIX + g.tx_name, g.chrom, start, end, 
                      g.gene_name, g.strand, seqlines))
        yield g, fa_record
    logging.info("Used %d/%d gene features" % (used,total))
    ref_fa.close()

def create_chimerascan_index(output_dir, 
                             genome_fasta_file, 
                             gene_feature_file,
                             bowtie_build_bin):
    # create output dir if it does not exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        logging.info("Created index directory: %s" % (output_dir))
    # copy reference fasta file to output dir and write 
    # gene features to destination directory
    index_fasta_file = os.path.join(output_dir, ALIGN_INDEX + ".fa")
    dst_gene_feature_file = os.path.join(output_dir, GENE_FEATURE_FILE)
    msg = "Building transcriptome sequences and gene features"
    if (up_to_date(index_fasta_file, genome_fasta_file) and
        up_to_date(index_fasta_file, gene_feature_file) and
        up_to_date(dst_gene_feature_file, gene_feature_file)):
        logging.info("[SKIPPED] %s" % (msg))
    else:
        logging.info(msg)
        # open the genome fasta file to check for an index
        logging.debug("Checking reference genome sequence file")
        fh = pysam.Fastafile(genome_fasta_file)
        fh.close()
        # write sequences from gene feature file
        logging.info("Adding transcript sequences and gene features to index")
        fasta_fh = open(index_fasta_file, "w")
        gene_fh = open(dst_gene_feature_file, "w")
        for g, fa_record in genepred_to_fasta(gene_feature_file, genome_fasta_file):
            print >>gene_fh, str(g)
            print >>fasta_fh, fa_record
        gene_fh.close()
        fasta_fh.close()
        # remove old fasta index
        if os.path.exists(index_fasta_file + ".fai"):
            os.remove(index_fasta_file + ".fai")
        # index the combined fasta file
        logging.info("Indexing the FASTA file")
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
    logging.info("Chimerascan index created successfully")
    return JOB_SUCCESS


def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = OptionParser("usage: %prog [options] <reference_genome.fa> "
                          "<genepred_genes.txt> <index_output_dir>")
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
                                       gene_feature_file, bowtie_build_bin)
    sys.exit(retcode)

if __name__ == '__main__':
    main()