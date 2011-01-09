'''
Created on Jan 5, 2011

@author: mkiyer
'''
import logging
import os
import sys
from optparse import OptionParser

import lib.config as config
from lib.config import JOB_SUCCESS, JOB_ERROR
from lib.base import check_executable, get_read_length
from lib.align import align

def check_command_line_args(options, args, parser):
    # check command line arguments
    if len(args) < 3:
        parser.error("Incorrect number of command line arguments")
    fastq_files = args[0:2]
    output_dir = args[2]
    # check that input fastq files exist
    read_lengths = []
    for mate,fastq_file in enumerate(fastq_files):
        if not os.path.isfile(args[0]):
            parser.error("mate '%d' fastq file '%s' is not valid" % 
                         (mate, fastq_file))
        logging.debug("Checking read length for file %s" % 
                      (fastq_file))
        read_lengths.append(get_read_length(fastq_file))
        logging.debug("Read length for file %s: %d" % 
                      (fastq_file, read_lengths[-1]))
    # check that mate read lengths are equal
    if len(set(read_lengths)) > 1:
        parser.error("read lengths mate1=%d and mate2=%d are unequal" % 
                     (read_lengths[0], read_lengths[1]))
    # check that seed length < read length
    if any(options.seed_length > rlen for rlen in read_lengths):
        parser.error("seed length %d cannot be longer than read length" % 
                     (options.seed_length))
    # check that output dir is not a regular file
    if os.path.exists(output_dir) and (not os.path.isdir(output_dir)):
        parser.error("Output directory name '%s' exists and is not a valid directory" % 
                     (output_dir))
    # check that samtools binary exists
    if check_executable(options.samtools_bin):
        logging.debug("Checking for 'samtools' binary... found")
    else:
        parser.error("samtools binary not found or not executable")
    # check that bowtie-build program exists
    if check_executable(options.bowtie_build_bin):
        logging.debug("Checking for 'bowtie-build' binary... found")
    else:
        parser.error("bowtie-build binary not found or not executable")
    # check that bowtie program exists
    if check_executable(options.bowtie_bin):
        logging.debug("Checking for 'bowtie' binary... found")
    else:
        parser.error("bowtie binary not found or not executable")
    # check that alignment index exists
    if os.path.isdir(options.index_dir):
        logging.debug("Checking for chimerascan index directory... found")
    else:
        parser.error("chimerascan alignment index directory '%s' not valid" % 
                     (options.index_dir))
    # check that alignment index file exists
    align_index_file = os.path.join(options.index_dir, config.BOWTIE_INDEX_FILE)
    if os.path.isfile(align_index_file):
        logging.debug("Checking for bowtie index file... found")
    else:
        parser.error("chimerascan bowtie index file '%s' invalid" % (align_index_file))
    # check for sufficient processors
    processors_needed = (config.BASE_PROCESSORS + 2*options.bowtie_threads)    
    if processors_needed > options.num_processors:
        parser.error("Not enough processor cores (approx. %d needed, %d available)" % 
                     (processors_needed, options.num_processors)) 

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = OptionParser("usage: %prog [options] <mate1.fq> <mate2.fq> <output_dir>")
    parser.add_option("-p", "--processors", dest="num_processors", 
                      type="int", default=config.MIN_PROCESSORS)
    parser.add_option("--samtools-bin", dest="samtools_bin", 
                      default="samtools", help="Path to 'samtools' program")
    parser.add_option("--bowtie-build-bin", dest="bowtie_build_bin", 
                      default="bowtie-build", 
                      help="Path to 'bowtie-build' program")
    parser.add_option("--bowtie-bin", dest="bowtie_bin", default="bowtie", 
                      help="Path to 'bowtie' program")
    parser.add_option("--bowtie-threads", dest="bowtie_threads", type="int",
                      default=1, help="Bowtie threads per mate")
    parser.add_option("--index", dest="index_dir",
                      help="Path to chimerascan index directory")
    parser.add_option("--multihits", type="int", dest="multihits", 
                      default=40)
    parser.add_option("--mismatches", type="int", dest="mismatches", 
                      default=2)
    parser.add_option("--seed-length", type="int", dest="seed_length", 
                      default=25)
    parser.add_option("--quals", dest="fastq_format")
    parser.add_option("--max-fragment-length", type="int", 
                      dest="max_fragment_length", default=600,
                      help="Largest expected fragment length (reads less"
                      " than this fragment length are assumed to be "
                      " genomically contiguous")    
    options, args = parser.parse_args()
    check_command_line_args(options, args, parser)
    # extract command line arguments
    fastq_files = args[0:2]
    output_dir = args[2]
    # create output dir if it does not exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        logging.info("Created index directory: %s" % (output_dir))    
    align_output_file = os.path.join(output_dir, config.DISCORDANT_READS_FILE)
    align_expr_file = os.path.join(output_dir, config.EXPRESSION_FILE)
    bowtie_index = os.path.join(options.index_dir, config.ALIGN_INDEX)
    gene_feature_file = os.path.join(options.index_dir, config.GENE_FEATURE_FILE)
    #
    # Discordant reads alignment step
    #
    try:
        align(align_output_file, align_expr_file, fastq_files, 
              options.seed_length,
              options.fastq_format, 
              options.multihits, 
              options.mismatches,
              options.bowtie_threads, 
              options.bowtie_bin, 
              bowtie_index, 
              gene_feature_file, 
              config.GENE_REF_PREFIX,
              options.max_fragment_length)
        retcode = JOB_SUCCESS
    except Exception:
        retcode = JOB_ERROR
    sys.exit(retcode)

if __name__ == '__main__':
    main()
