'''
Created on Jan 5, 2011

@author: mkiyer
'''
import logging
import os
import sys
import subprocess
from optparse import OptionParser

import pysam

# local imports
import lib.config as config
from lib.config import JOB_SUCCESS, JOB_ERROR
from lib.base import check_executable, get_read_length, parse_library_type
from lib.align import align
from lib.merge_read_pairs import merge_read_pairs
from lib.find_discordant_reads import discordant_reads_to_chimeras
from lib.nominate_chimeras import nominate_chimeras
from lib.bedpe_to_fasta import bedpe_to_junction_fasta
from lib.merge_spanning_alignments import merge_spanning_alignments
from lib.sort_discordant_reads import sort_discordant_reads

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
    if any(options.segment_length > rlen for rlen in read_lengths):
        parser.error("seed length %d cannot be longer than read length" % 
                     (options.segment_length))
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
    if options.num_processors < config.BASE_PROCESSORS:
        logging.warning("Please specify >=2 processes using '-p' to allow program to run efficiently")

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = OptionParser("usage: %prog [options] <mate1.fq> <mate2.fq> <output_dir>")
    parser.add_option("-p", "--processors", dest="num_processors", 
                      type="int", default=config.BASE_PROCESSORS)
    parser.add_option("--index", dest="index_dir",
                      help="Path to chimerascan index directory")
    parser.add_option("--samtools-bin", dest="samtools_bin", 
                      default="samtools", help="Path to 'samtools' program")
    parser.add_option("--bowtie-build-bin", dest="bowtie_build_bin", 
                      default="bowtie-build", 
                      help="Path to 'bowtie-build' program")
    parser.add_option("--bowtie-bin", dest="bowtie_bin", default="bowtie", 
                      help="Path to 'bowtie' program")
    parser.add_option("--bowtie-mode-v", action="store_true", 
                      dest="bowtie_mode_v", default=False,
                      help="Run bowtie with -v to ignore quality scores")
    parser.add_option("--multihits", type="int", dest="multihits", 
                      default=40)
    parser.add_option("--mismatches", type="int", dest="mismatches", 
                      default=2)
    parser.add_option("--segment-length", type="int", dest="segment_length", 
                      default=25)
    parser.add_option("--trim5", type="int", dest="trim5", 
                      default=0)
    parser.add_option("--trim3", type="int", dest="trim3", 
                      default=0)    
    parser.add_option("--quals", dest="fastq_format")
    parser.add_option("--min-fragment-length", type="int", 
                      dest="min_fragment_length", default=50,
                      help="Smallest expected fragment length")
    parser.add_option("--max-fragment-length", type="int", 
                      dest="max_fragment_length", default=600,
                      help="Largest expected fragment length (reads less"
                      " than this fragment length are assumed to be "
                      " genomically contiguous")
    parser.add_option('--library', dest="library_type", default="fr")    
    options, args = parser.parse_args()
    check_command_line_args(options, args, parser)
    # extract command line arguments
    fastq_files = args[0:2]
    output_dir = args[2]
    # create output dir if it does not exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        logging.info("Created index directory: %s" % (output_dir))    
    library_type = parse_library_type(options.library_type)    
    gene_feature_file = os.path.join(options.index_dir, config.GENE_FEATURE_FILE)
    try:
        #
        # Discordant reads alignment step
        #
        align_output_file = os.path.join(output_dir, config.ALIGNED_READS_FILE)
        bowtie_index = os.path.join(options.index_dir, config.ALIGN_INDEX)        
        align(fastq_files, options.fastq_format, bowtie_index,
              align_output_file, options.bowtie_bin, 
              options.num_processors, options.segment_length,
              options.trim5, options.trim3, options.multihits,
              options.mismatches, options.bowtie_mode)        
        #
        # Merge paired-end reads step
        #
        paired_bam_file = os.path.join(output_dir, config.PAIRED_ALIGNED_READS_FILE)        
        bamfh = pysam.Samfile(align_output_file, "rb")
        paired_bamfh = pysam.Samfile(paired_bam_file, "wb", template=bamfh)
        merge_read_pairs(bamfh, paired_bamfh, 
                         options.min_fragment_length,
                         options.max_fragment_length,
                         library_type)
        paired_bamfh.close() 
        bamfh.close()
        #
        # Find discordant reads step
        #
        discordant_bedpe_file = os.path.join(output_dir, config.DISCORDANT_BEDPE_FILE)
        unmapped_fasta_file = os.path.join(output_dir, config.UNMAPPED_FASTA_FILE)
        # TODO: add contam refs
        bamfh = pysam.Samfile(paired_bam_file, "rb")
        discordant_reads_to_chimeras(bamfh, discordant_bedpe_file, gene_feature_file,
                                     options.max_fragment_length, library_type,
                                     contam_refs=None,
                                     unmapped_fasta_file=unmapped_fasta_file)
        bamfh.close()
        #
        # Sort discordant reads
        #
        sorted_discordant_bedpe_file = os.path.join(output_dir, config.SORTED_DISCORDANT_BEDPE_FILE)
        sort_discordant_reads(discordant_bedpe_file, sorted_discordant_bedpe_file)        
        #
        # Nominate chimeras step
        #
        encompassing_bedpe_file = os.path.join(output_dir, config.ENCOMPASSING_CHIMERA_BEDPE_FILE)        
        infh = open(sorted_discordant_bedpe_file, "r")
        outfh = open(encompassing_bedpe_file, "w")                
        # TODO: add contam refs
        nominate_chimeras(infh, outfh, gene_feature_file,
                          contam_refs=None, 
                          trim=config.EXON_JUNCTION_TRIM_BP)
        outfh.close()
        infh.close()        
        #
        # Extract junction sequences from chimeras file
        #        
        ref_fasta_file = os.path.join(options.index_dir, config.ALIGN_INDEX + ".fa")
        read_length = get_read_length(fastq_files[0])
        junc_fasta_file = os.path.join(options.output_dir, config.JUNC_REF_FASTA_FILE)
        junc_map_file = os.path.join(options.output_dir, config.JUNC_REF_MAP_FILE)        
        bedpe_to_junction_fasta(encompassing_bedpe_file, ref_fasta_file,                                
                                read_length, open(junc_fasta_file, "w"),
                                open(junc_map_file, "w"))
        #
        # Build a bowtie index to align and detect spanning reads
        #
        bowtie_spanning_index = os.path.join(options.output_dir, config.JUNC_BOWTIE_INDEX)
        args = [options.bowtie_build_bin, junc_fasta_file, bowtie_spanning_index]
        subprocess.call(args)
        #
        # Align unmapped reads across putative junctions
        #
        junc_bam_file = os.path.join(output_dir, config.JUNC_READS_BAM_FILE)
        bowtie_index = os.path.join(options.index_dir, config.ALIGN_INDEX)        
        align([unmapped_fasta_file], "phred33-quals", bowtie_spanning_index,
              junc_bam_file, options.bowtie_bin, 
              options.num_processors, options.segment_length,
              options.trim5, options.trim3, options.multihits,
              options.mismatches, options.bowtie_mode)
        #
        # Merge spanning and encompassing read information
        #
        chimera_bedpe_file = os.path.join(output_dir, config.CHIMERA_BEDPE_FILE)        
        merge_spanning_alignments(junc_bam_file, junc_map_file, chimera_bedpe_file,
                                  read_length, anchor_min=0, anchor_max=0,
                                  anchor_mismatches=0)
        retcode = JOB_SUCCESS
    except Exception:
        retcode = JOB_ERROR
    sys.exit(retcode)

if __name__ == '__main__':
    main()
