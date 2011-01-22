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
from lib.align_segments import align
from lib.align_full import align_pe_full
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

def up_to_date(outfile, infile):
    if not os.path.exists(infile):
        return False
    if not os.path.exists(outfile):
        return False
    if os.path.getsize(outfile) == 0:
        return False    
    return os.path.getmtime(outfile) >= os.path.getmtime(infile)

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
        logging.info("Created output directory: %s" % (output_dir))    
    # gather and parse run parameters
    library_type = parse_library_type(options.library_type)    
    gene_feature_file = os.path.join(options.index_dir, config.GENE_FEATURE_FILE)
    bowtie_mode = "-v" if options.bowtie_mode_v else "-n"
    bowtie_index = os.path.join(options.index_dir, config.ALIGN_INDEX)
    read_length = get_read_length(fastq_files[0])
    # minimum fragment length cannot be smaller than the trimmed read length
    trimmed_read_length = read_length - options.trim5 - options.trim3
    min_fragment_length = max(options.min_fragment_length, 
                              trimmed_read_length)
    #
    # Initial Bowtie alignment step
    #
    # align in paired-end mode, trying to resolve as many reads as possible
    # this effectively rules out the vast majority of reads as candidate
    # fusions
    unaligned_fastq_param = os.path.join(output_dir, config.UNALIGNED_FASTQ_PARAM)
    maxmultimap_fastq_param = os.path.join(output_dir, config.MAXMULTIMAP_FASTQ_PARAM)
    aligned_bam_file = os.path.join(output_dir, config.ALIGNED_READS_BAM_FILE)
    if all(up_to_date(aligned_bam_file, fq) for fq in fastq_files):
        logging.info("[SKIPPED] Alignment results exist")
    else:    
        logging.info("Aligning full-length reads in paired-end mode")
        retcode = align_pe_full(fastq_files, 
                                bowtie_index,
                                aligned_bam_file, 
                                unaligned_fastq_param,
                                maxmultimap_fastq_param,
                                min_fragment_length=min_fragment_length,
                                max_fragment_length=options.max_fragment_length,
                                trim5=options.trim5,
                                trim3=options.trim3,
                                library_type=options.library_type,
                                num_processors=options.num_processors,
                                fastq_format=options.fastq_format,
                                multihits=options.multihits,
                                mismatches=options.mismatches,
                                bowtie_bin=options.bowtie_bin,
                                bowtie_mode=bowtie_mode)
        if retcode != 0:
            logging.error("Bowtie failed with error code %d" % (retcode))    
            sys.exit(retcode)
    #
    # Discordant reads alignment step
    #
    discordant_bam_file = os.path.join(output_dir, config.DISCORDANT_BAM_FILE)
    unaligned_fastq_files = [os.path.join(output_dir, fq) for fq in config.UNALIGNED_FASTQ_FILES]    
    if all(up_to_date(discordant_bam_file, fq) for fq in fastq_files):
        logging.info("[SKIPPED] Discordant alignment results exist")
    else:
        logging.info("Aligning initially unmapped reads in single read mode")
        align(unaligned_fastq_files, options.fastq_format, bowtie_index,
              discordant_bam_file, 
              bowtie_bin=options.bowtie_bin,
              num_processors=options.num_processors, 
              segment_length=options.segment_length,
              segment_trim=False,
              trim5=options.trim5, 
              trim3=options.trim3, 
              multihits=options.multihits,
              mismatches=options.mismatches, 
              bowtie_mode=bowtie_mode)
    #
    # Merge paired-end reads step
    #
    paired_bam_file = os.path.join(output_dir, config.DISCORDANT_PAIRED_BAM_FILE)
    if up_to_date(paired_bam_file, discordant_bam_file):
        logging.info("[SKIPPED] Read pairing results exist")
    else:
        logging.info("Pairing aligned reads")
        bamfh = pysam.Samfile(discordant_bam_file, "rb")
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
    spanning_fastq_file = os.path.join(output_dir, config.SPANNING_FASTQ_FILE)
    if (up_to_date(discordant_bedpe_file, paired_bam_file) and
        up_to_date(spanning_fastq_file, paired_bam_file)):    
        logging.info("[SKIPPED] Discordant BEDPE file exists")
    else:
        logging.info("Nominating discordant reads")
        # TODO: add contam refs
        discordant_reads_to_chimeras(paired_bam_file, discordant_bedpe_file, gene_feature_file,
                                     options.max_fragment_length, library_type,
                                     contam_refs=None,
                                     unmapped_fastq_file=spanning_fastq_file)
    #
    # Sort discordant reads
    #
    sorted_discordant_bedpe_file = os.path.join(output_dir, config.SORTED_DISCORDANT_BEDPE_FILE)
    if (up_to_date(sorted_discordant_bedpe_file, discordant_bedpe_file)):
        logging.info("[SKIPPED] Sorted discordant BEDPE file exists")
    else:        
        logging.info("Sorting discordant reads")
        sort_discordant_reads(discordant_bedpe_file, sorted_discordant_bedpe_file)        
    #
    # Nominate chimeras step
    #
    encompassing_bedpe_file = os.path.join(output_dir, config.ENCOMPASSING_CHIMERA_BEDPE_FILE)        
    if (up_to_date(encompassing_bedpe_file, sorted_discordant_bedpe_file)):
        logging.info("[SKIPPED] Encompassing chimeras BEDPE file exists")
    else:        
        logging.info("Aggregating discordant reads across genes")
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
    junc_fasta_file = os.path.join(output_dir, config.JUNC_REF_FASTA_FILE)
    junc_map_file = os.path.join(output_dir, config.JUNC_REF_MAP_FILE)
    if (up_to_date(junc_fasta_file, encompassing_bedpe_file) and
        up_to_date(junc_map_file, encompassing_bedpe_file)):        
        logging.info("[SKIPPED] Chimeric junction files exist")
    else:        
        logging.info("Extracting chimeric junction sequences")
        bedpe_to_junction_fasta(encompassing_bedpe_file, ref_fasta_file,                                
                                read_length, open(junc_fasta_file, "w"),
                                open(junc_map_file, "w"))
    #
    # Build a bowtie index to align and detect spanning reads
    #
    bowtie_spanning_index = os.path.join(output_dir, config.JUNC_BOWTIE_INDEX)
    if (up_to_date(bowtie_spanning_index, junc_fasta_file)):
        logging.info("[SKIPPED] Bowtie junction index exists")
    else:        
        logging.info("Building bowtie index for junction-spanning reads")
        args = [options.bowtie_build_bin, junc_fasta_file, bowtie_spanning_index]
        subprocess.call(args)
    #
    # Align unmapped reads across putative junctions
    #
    junc_bam_file = os.path.join(output_dir, config.JUNC_READS_BAM_FILE)
    if (up_to_date(junc_bam_file, bowtie_spanning_index) and
        up_to_date(junc_bam_file, spanning_fastq_file)):
        logging.info("[SKIPPED] Spanning read alignment exists")
    else:            
        logging.info("Aligning junction-spanning reads")
        align([spanning_fastq_file], 
              "phred33-quals", 
              bowtie_spanning_index,
              junc_bam_file, 
              bowtie_bin=options.bowtie_bin, 
              num_processors=options.num_processors, 
              segment_length=options.segment_length,
              segment_trim=False,
              trim5=options.trim5, 
              trim3=options.trim3, 
              multihits=options.multihits,
              mismatches=options.mismatches, 
              bowtie_mode=bowtie_mode)
    #
    # Merge spanning and encompassing read information
    #
    chimera_bedpe_file = os.path.join(output_dir, config.CHIMERA_BEDPE_FILE)
    if (up_to_date(chimera_bedpe_file, junc_bam_file) and
        up_to_date(chimera_bedpe_file, junc_map_file)):
        logging.info("[SKIPPED] Chimera BEDPE file exists")
    else:
        logging.info("Merging spanning and encompassing read alignments")
        merge_spanning_alignments(junc_bam_file, junc_map_file, chimera_bedpe_file,
                                  read_length, anchor_min=0, anchor_max=0,
                                  anchor_mismatches=0)
    retcode = JOB_SUCCESS
    sys.exit(retcode)

if __name__ == '__main__':
    main()
