'''
Created on Jan 5, 2011

@author: mkiyer
'''
__author__ = "Matthew Iyer"
__copyright__ = "Copyright 2011, chimerascan project"
__credits__ = ["Matthew Iyer", "Christopher Maher"]
__license__ = "GPL"
__version__ = "0.3.0"
__maintainer__ = "Matthew Iyer"
__email__ = "mkiyer@med.umich.edu"
__status__ = "beta"

import logging
import os
import subprocess
import sys
from optparse import OptionParser, OptionGroup
import xml.etree.ElementTree as etree

# check for python version 2.6.0 or greater
if sys.version_info < (2,6,0):
    sys.stderr.write("You need python 2.6 or later to run chimerascan\n")
    sys.exit(1)

# local imports
import pysam
import lib.config as config
from lib.config import JOB_SUCCESS, JOB_ERROR
from lib.base import check_executable, get_read_length, parse_library_type, parse_bool

from pipeline.align_full import align_pe_full, align_sr_full
from pipeline.align_segments import align, determine_read_segments
from pipeline.merge_read_pairs import merge_read_pairs
from pipeline.find_discordant_reads import find_discordant_reads
from pipeline.extend_sequences import extend_sequences
from pipeline.sort_discordant_reads import sort_discordant_reads
from pipeline.nominate_chimeras import nominate_chimeras
from pipeline.filter_encompassing_chimeras import filter_encompassing_chimeras
from pipeline.nominate_spanning_reads import nominate_spanning_reads
from pipeline.bedpe_to_fasta import bedpe_to_junction_fasta
from pipeline.merge_spanning_alignments import merge_spanning_alignments
from pipeline.profile_insert_size import profile_isize_stats
from pipeline.filter_spanning_chimeras import filter_spanning_chimeras


DEFAULT_NUM_PROCESSORS = config.BASE_PROCESSORS
DEFAULT_KEEP_TMP = False
DEFAULT_BOWTIE_BUILD_BIN = "bowtie-build"
DEFAULT_BOWTIE_BIN = "bowtie"
DEFAULT_BOWTIE_MODE_V = False
DEFAULT_MULTIHITS = 100
DEFAULT_MISMATCHES = 2
DEFAULT_SEGMENT_LENGTH = 25
DEFAULT_TRIM5 = 0
DEFAULT_TRIM3 = 0
DEFAULT_MIN_FRAG_LENGTH = 0
DEFAULT_MAX_FRAG_LENGTH = 1000 
DEFAULT_MAX_INDEL_SIZE = 100
DEFAULT_FASTQ_FORMAT = "sanger"
DEFAULT_LIBRARY_TYPE = "fr"

DEFAULT_FILTER_MAX_MULTIMAPS = 2
DEFAULT_FILTER_MULTIMAP_RATIO = 0.10
DEFAULT_FILTER_STRAND_PVALUE = 0.01
DEFAULT_FILTER_ISIZE_STDEVS = 3

DEFAULT_ANCHOR_MIN = 0
DEFAULT_ANCHOR_MAX = 5
DEFAULT_ANCHOR_MISMATCHES = 5

translate_quals = {'solexa': 'solexa-quals',
                   'illumina': 'solexa1.3-quals',
                   'sanger': 'phred33-quals'}
quals_help_string = ",".join(translate_quals.keys())

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

class RunConfig(object):
    def __init__(self):
        self.output_dir = None
        self.fastq_files = None
        self.index_dir = None
        self.fastq_format = None 
        self.num_processors = None
        self.keep_tmp = None
        self.bowtie_build_bin = None
        self.bowtie_bin = None
        self.bowtie_mode_v = None
        self.multihits = None
        self.mismatches = None
        self.segment_length = None
        self.trim5 = None
        self.trim3 = None
        self.min_fragment_length = None
        self.max_fragment_length = None
        self.max_indel_size = None
        self.library_type = None
        # filtering params
        self.filter_max_multimaps = None
        self.filter_multimap_ratio = None
        self.filter_isize_stdevs = None
        self.filter_strand_pval = None
        # spanning read constraints
        self.anchor_min = None
        self.anchor_max = None 
        self.anchor_mismatches = None
        

    def from_xml(self, xmlfile):
        tree = etree.parse(xmlfile)        
        root = tree.getroot()
        # required arguments        
        self.output_dir = root.findtext('output_dir')
        fastq_files = {}
        for mate_elem in root.findall("fastq_files/*"):
            mate = int(mate_elem.get("mate"))
            fastq_files[mate] = mate_elem.text
        self.fastq_files = [fastq_files[mate] for mate in xrange(len(fastq_files))]
        self.index_dir = root.findtext('index')        
        # optional arguments
        quals = root.findtext("quals")
        if quals not in translate_quals:
            logging.error("Quality score option %s unknown, using default %s" % 
                          (quals, DEFAULT_FASTQ_FORMAT))
        self.fastq_format = translate_quals.get(quals, DEFAULT_FASTQ_FORMAT)
        self.num_processors = int(root.findtext("num_processors"), DEFAULT_NUM_PROCESSORS)
        self.keep_tmp = parse_bool(root.findtext("keep_tmp", str(DEFAULT_KEEP_TMP)))        
        self.bowtie_build_bin = root.findtext("bowtie_build_bin", DEFAULT_BOWTIE_BUILD_BIN)
        self.bowtie_bin = root.findtext("bowtie_bin", DEFAULT_BOWTIE_BIN)
        self.bowtie_mode_v = parse_bool(root.findtext("bowtie_mode_v", str(DEFAULT_BOWTIE_MODE_V)))
        self.multihits = int(root.findtext("multihits", DEFAULT_MULTIHITS))
        self.mismatches = int(root.findtext("mismatches", DEFAULT_MISMATCHES))
        self.segment_length = int(root.findtext("segment_length", DEFAULT_SEGMENT_LENGTH))
        self.trim5 = root.findtext("trim5", DEFAULT_TRIM5)
        self.trim3 = root.findtext("trim3", DEFAULT_TRIM3)
        self.min_fragment_length = root.findtext("min_fragment_length", DEFAULT_MIN_FRAG_LENGTH)
        self.max_fragment_length = root.findtext("max_fragment_length", DEFAULT_MAX_FRAG_LENGTH)
        self.max_indel_size = root.findtext("max_indel_size", DEFAULT_MAX_INDEL_SIZE)
        self.library_type = root.findtext("library_type", DEFAULT_LIBRARY_TYPE)
        # filtering params
        self.filter_max_multimaps = root.findtext("filter_max_multimaps", DEFAULT_FILTER_MAX_MULTIMAPS)
        self.filter_multimap_ratio = root.findtext("filter_multimap_ratio", DEFAULT_FILTER_MULTIMAP_RATIO)
        self.filter_isize_stdevs = root.findtext("filter_isize_stdevs", DEFAULT_FILTER_ISIZE_STDEVS)
        self.filter_strand_pval = root.findtext("filter_strand_pval", DEFAULT_FILTER_STRAND_PVALUE)
        # spanning read constraints
        self.anchor_min = root.findtext("anchor_min", DEFAULT_ANCHOR_MIN)
        self.anchor_max = root.findtext("anchor_min", DEFAULT_ANCHOR_MAX)
        self.anchor_mismatches = root.findtext("anchor_mismatches", DEFAULT_ANCHOR_MISMATCHES)

    def from_args(self, args):
        parser = OptionParser(usage="%prog [options] [--config <config_file> "
                              " | <mate1.fq> <mate2.fq> <output_dir>]",
                              version="%s" % __version__)
        # standard options
        parser.add_option("--index", dest="index_dir", default=None,
                          help="Path to chimerascan index directory")
        parser.add_option("--config", dest="config_file",
                          help="Path to configuration XML file") 
        parser.add_option("-p", "--processors", dest="num_processors", 
                          type="int", default=DEFAULT_NUM_PROCESSORS,
                          help="Number of processor cores to allocate to "
                          "chimerascan [default=%default]")
        parser.add_option("--keep-tmp", dest="keep_tmp", action="store_true", 
                          default=DEFAULT_KEEP_TMP,
                          help="Do not delete intermediate files from run") 
        parser.add_option("--bowtie-build-bin", dest="bowtie_build_bin", 
                          default=DEFAULT_BOWTIE_BUILD_BIN, 
                          help="Path to 'bowtie-build' program")
        parser.add_option("--bowtie-bin", dest="bowtie_bin", default=DEFAULT_BOWTIE_BIN, 
                          help="Path to 'bowtie' program")
        parser.add_option("--bowtie-mode-v", action="store_true", 
                          dest="bowtie_mode_v", default=DEFAULT_BOWTIE_MODE_V,
                          help="Run bowtie with -v to ignore quality scores")
        parser.add_option("--multihits", type="int", dest="multihits", 
                          default=DEFAULT_MULTIHITS, metavar="X",
                          help="Ignore reads that map to more than X "
                          "locations [default=%default]")
        parser.add_option("--mismatches", type="int", dest="mismatches",
                          default=DEFAULT_MISMATCHES, metavar="X",
                          help="Aligned reads must have <= X mismatches "
                          "[default=%default]")
        parser.add_option("--segment-length", type="int", dest="segment_length", 
                          default=DEFAULT_SEGMENT_LENGTH,
                          help="Size of read segments during discordant " 
                          "alignment phase [default=%default]")
        parser.add_option("--trim5", type="int", dest="trim5", 
                          default=DEFAULT_TRIM5, metavar="N",
                          help="Trim N bases from 5' end of read")
        parser.add_option("--trim3", type="int", dest="trim3", 
                          default=DEFAULT_TRIM3, metavar="N",
                          help="Trim N bases from 3' end of read")
        parser.add_option("--quals", dest="fastq_format", 
                          default=DEFAULT_FASTQ_FORMAT, metavar="FMT",
                          help="Choose from %s [default=%s]" % 
                          (quals_help_string, DEFAULT_FASTQ_FORMAT))
        parser.add_option("--min-fragment-length", type="int", 
                          dest="min_fragment_length", 
                          default=DEFAULT_MIN_FRAG_LENGTH,
                          help="Smallest expected fragment length "
                          "[default=%default]")
        parser.add_option("--max-fragment-length", type="int", 
                          dest="max_fragment_length", 
                          default=DEFAULT_MAX_FRAG_LENGTH,
                          help="Largest expected fragment length (reads less"
                          " than this fragment length are assumed to be "
                          " genomically contiguous) [default=%default]")
        parser.add_option("--max-indel-size", type="int", 
                          dest="max_indel_size", 
                          default=DEFAULT_MAX_INDEL_SIZE,
                          help="Tolerate indels less than N bp "
                          "[default=%default]", metavar="N")
        parser.add_option('--library', dest="library_type", 
                          default=DEFAULT_LIBRARY_TYPE,
                          help="Library type ('fr', 'rf') [default=%default]")
        parser.add_option("--anchor-min", type="int", dest="anchor_min", 
                          default=DEFAULT_ANCHOR_MIN,
                          help="Minimum junction overlap required to call "
                          "spanning reads")
        parser.add_option("--anchor-max", type="int", dest="anchor_max", 
                          default=DEFAULT_ANCHOR_MAX,
                          help="Junction overlap below which to enforce "
                          "mismatch checks")
        parser.add_option("--anchor-mismatches", type="int", dest="anchor_mismatches", 
                          default=DEFAULT_ANCHOR_MISMATCHES,
                          help="Number of mismatches allowed within anchor "
                          "region")
        filter_group = OptionGroup(parser, "Filtering options",
                                   "Adjust these options to change "
                                   "filtering behavior") 
        filter_group.add_option("--filter-multimaps", type="int",
                                dest="filter_max_multimaps",
                                default=DEFAULT_FILTER_MAX_MULTIMAPS,
                                help="Filter chimeras that lack a read "
                                " with <= HITS alternative hits in both "
                                " pairs [default=%default]")
        filter_group.add_option("--filter-multimap-ratio", type="float",
                                default=DEFAULT_FILTER_MULTIMAP_RATIO,
                                dest="filter_multimap_ratio", metavar="RATIO",
                                help="Filter chimeras with a weighted coverage "
                                "versus total reads ratio <= RATIO "
                                "[default=%default]")
        filter_group.add_option("--filter-isize-stdevs", type="int",
                                default=DEFAULT_FILTER_ISIZE_STDEVS,
                                dest="filter_isize_stdevs", metavar="N",
                                help="Filter chimeras where putative insert "
                                "size is >N standard deviations from the "
                                "mean [default=%default]")
        filter_group.add_option("--filter-strand-pvalue", type="float",
                                default=DEFAULT_FILTER_STRAND_PVALUE,
                                dest="filter_strand_pval", metavar="p",
                                help="p-value to reject chimera based on "
                                " binomial test that balance of +/- strand "
                                " encompassing reads should be 50/50 "
                                "[default=%default]")
        parser.add_option_group(filter_group)        
        options, args = parser.parse_args(args=args)
        # parse config file options/args
        if options.config_file is not None:
            self.from_xml(options.config_file)
        # check command line arguments
        if self.fastq_files is None:
            if len(args) < 2:
                parser.error("fastq files not specified in config file or command line")        
            self.fastq_files = args[0:2]
        if self.output_dir is None:
            if len(args) < 3:
                parser.error("output dir not specified in config file or command line")   
            self.output_dir = args[2]
        if self.index_dir is None:
            if options.index_dir is None:
                parser.error("index dir not specified in config file or command line")                   
            self.index_dir = options.index_dir
        # optional arguments
        if self.fastq_format is None:
            self.fastq_format = translate_quals.get(options.fastq_format, DEFAULT_FASTQ_FORMAT)
        if self.num_processors is None:
            self.num_processors = int(options.num_processors)
        if self.keep_tmp is None:
            self.keep_tmp = options.keep_tmp
        if self.bowtie_build_bin is None:
            self.bowtie_build_bin = options.bowtie_build_bin
        if self.bowtie_bin is None:
            self.bowtie_bin = options.bowtie_bin
        if self.bowtie_mode_v is None:
            self.bowtie_mode_v = options.bowtie_mode_v
        if self.multihits is None:
            self.multihits = options.multihits
        if self.mismatches is None:
            self.mismatches = options.mismatches
        if self.segment_length is None:
            self.segment_length = options.segment_length
        if self.trim5 is None:
            self.trim5 = options.trim5
        if self.trim3 is None:
            self.trim3 = options.trim3
        if self.min_fragment_length is None:
            self.min_fragment_length = options.min_fragment_length
        if self.max_fragment_length is None:
            self.max_fragment_length = options.max_fragment_length
        if self.max_indel_size is None:
            self.max_indel_size = options.max_indel_size
        if self.library_type is None:
            self.library_type = options.library_type
        # set rest of options, overriding if attribute is undefined
        # or set to something other than the default 
        attrs = (("filter_max_multimaps", DEFAULT_FILTER_MAX_MULTIMAPS),
                 ("filter_multimap_ratio", DEFAULT_FILTER_MULTIMAP_RATIO),
                 ("filter_isize_stdevs", DEFAULT_FILTER_ISIZE_STDEVS),
                 ("filter_strand_pval", DEFAULT_FILTER_STRAND_PVALUE),
                 ("anchor_min", DEFAULT_ANCHOR_MIN),
                 ("anchor_max", DEFAULT_ANCHOR_MAX),
                 ("anchor_mismatches", DEFAULT_ANCHOR_MISMATCHES))
        for attr_name,default_val in attrs:
            if ((getattr(self, attr_name) is None) or
                (getattr(options, attr_name) != default_val)):
                setattr(self, attr_name, getattr(options, attr_name)) 

    def check_config(self):
        # check that input fastq files exist
        config_passed = True
        read_lengths = []
        for mate,fastq_file in enumerate(self.fastq_files):
            if not os.path.isfile(fastq_file):
                logging.error("mate '%d' fastq file '%s' is not valid" % 
                              (mate, fastq_file))
                config_passed = False
            read_lengths.append(get_read_length(fastq_file))
            logging.debug("Checking file %s" % (fastq_file))
            logging.debug("File %s read length=%d" % (fastq_file, read_lengths[-1]))
        # check that mate read lengths are equal
        if len(set(read_lengths)) > 1:
            logging.error("Unequal read lengths mate1=%d and mate2=%d" % 
                          (read_lengths[0], read_lengths[1]))
            config_passed = False
        # check that seed length < read length
        if any(self.segment_length > rlen for rlen in read_lengths):
            logging.error("seed length %d cannot be longer than read length" % 
                         (self.segment_length))
            config_passed = False
        # check that output dir is not a regular file
        if os.path.exists(self.output_dir) and (not os.path.isdir(self.output_dir)):
            logging.error("Output directory name '%s' exists and is not a valid directory" % 
                          (self.output_dir))
            config_passed = False
        if check_executable(self.bowtie_build_bin):
            logging.debug("Checking for 'bowtie-build' binary... found")
        else:
            logging.error("bowtie-build binary not found or not executable")
            config_passed = False
        # check that bowtie program exists
        if check_executable(self.bowtie_bin):
            logging.debug("Checking for 'bowtie' binary... found")
        else:
            logging.error("bowtie binary not found or not executable")
            config_passed = False
        # check that alignment index exists
        if os.path.isdir(self.index_dir):
            logging.debug("Checking for chimerascan index directory... found")
            # check that alignment index file exists
            align_index_file = os.path.join(self.index_dir, config.BOWTIE_INDEX_FILE)
            if os.path.isfile(align_index_file):
                logging.debug("Checking for bowtie index file... found")
            else:
                logging.error("chimerascan bowtie index file '%s' invalid" % (align_index_file))
                config_passed = False
        else:
            logging.error("chimerascan alignment index directory '%s' not valid" % 
                          (self.index_dir))
            config_passed = False
        # check for sufficient processors
        if self.num_processors < config.BASE_PROCESSORS:
            logging.warning("Please specify >=2 processes using '-p' to allow program to run efficiently")
        return config_passed

def up_to_date(outfile, infile):
    if not os.path.exists(infile):
        return False
    if not os.path.exists(outfile):
        return False
    if os.path.getsize(outfile) == 0:
        return False    
    return os.path.getmtime(outfile) >= os.path.getmtime(infile)

def run_chimerascan(runconfig):
    # normal run
    config_passed = runconfig.check_config()
    if not config_passed:
        logging.error("Invalid run configuration, aborting.")
        sys.exit(JOB_ERROR)
    # create output dir if it does not exist
    if not os.path.exists(runconfig.output_dir):
        os.makedirs(runconfig.output_dir)
        logging.info("Created output directory: %s" % (runconfig.output_dir))
    # create log dir if it does not exist
    log_dir = os.path.join(runconfig.output_dir, config.LOG_DIR)
    if not os.path.exists(log_dir):
        os.makedirs(log_dir)
        logging.debug("Created directory for log files: %s" % (log_dir))        
    # create tmp dir if it does not exist
    tmp_dir = os.path.join(runconfig.output_dir, config.TMP_DIR)
    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir)
        logging.debug("Created directory for tmp files: %s" % (tmp_dir))
    # gather and parse run parameters
    library_type = parse_library_type(runconfig.library_type)    
    gene_feature_file = os.path.join(runconfig.index_dir, config.GENE_FEATURE_FILE)
    bowtie_mode = "-v" if runconfig.bowtie_mode_v else "-n"
    bowtie_index = os.path.join(runconfig.index_dir, config.ALIGN_INDEX)
    original_read_length = get_read_length(runconfig.fastq_files[0])
    # minimum fragment length cannot be smaller than the trimmed read length
    trimmed_read_length = original_read_length - runconfig.trim5 - runconfig.trim3
    min_fragment_length = max(runconfig.min_fragment_length, 
                              trimmed_read_length)
    #
    # Initial Bowtie alignment step
    #
    # align in paired-end mode, trying to resolve as many reads as possible
    # this effectively rules out the vast majority of reads as candidate
    # fusions
    unaligned_fastq_param = os.path.join(tmp_dir, config.UNALIGNED_FASTQ_PARAM)
    maxmultimap_fastq_param = os.path.join(tmp_dir, config.MAXMULTIMAP_FASTQ_PARAM)
    aligned_bam_file = os.path.join(runconfig.output_dir, config.ALIGNED_READS_BAM_FILE)
    if all(up_to_date(aligned_bam_file, fq) for fq in runconfig.fastq_files):
        logging.info("[SKIPPED] Alignment results exist")
    else:    
        logging.info("Aligning full-length reads in paired-end mode")
        retcode = align_pe_full(runconfig.fastq_files, 
                                bowtie_index,
                                aligned_bam_file, 
                                unaligned_fastq_param,
                                maxmultimap_fastq_param,
                                min_fragment_length=min_fragment_length,
                                max_fragment_length=runconfig.max_fragment_length,
                                trim5=runconfig.trim5,
                                trim3=runconfig.trim3,
                                library_type=runconfig.library_type,
                                num_processors=runconfig.num_processors,
                                fastq_format=runconfig.fastq_format,
                                multihits=runconfig.multihits,
                                mismatches=runconfig.mismatches,
                                bowtie_bin=runconfig.bowtie_bin,
                                bowtie_mode=bowtie_mode)
        if retcode != 0:
            logging.error("Bowtie failed with error code %d" % (retcode))    
            sys.exit(retcode)
    #
    # Get insert size distribution
    #
    isize_stats_file = os.path.join(runconfig.output_dir, config.ISIZE_STATS_FILE)
    if up_to_date(isize_stats_file, aligned_bam_file):
        logging.info("[SKIPPED] Profiling insert size distribution")
        f = open(isize_stats_file, "r")
        isize_stats = f.next().strip().split('\t')
        f.close()
    else:
        logging.info("Profiling insert size distribution")
        max_isize_samples = config.ISIZE_MAX_SAMPLES
        bamfh = pysam.Samfile(aligned_bam_file, "rb")
        isize_stats = profile_isize_stats(bamfh,
                                          min_isize=min_fragment_length,
                                          max_isize=runconfig.max_fragment_length,
                                          max_samples=max_isize_samples)
    # unpack insert size statistics tuple for use in downstream stages
    isize_mean, isize_median, isize_mode, isize_std = isize_stats    
    #
    # Discordant reads alignment step
    #
    discordant_bam_file = os.path.join(tmp_dir, config.DISCORDANT_BAM_FILE)
    unaligned_fastq_files = [os.path.join(tmp_dir, fq) for fq in config.UNALIGNED_FASTQ_FILES]
    # get the segments used in discordant alignment to know the effective
    # read length used to align.  we used this to set the 'padding' during
    # spanning read discovery
    segments = determine_read_segments(original_read_length, 
                                       segment_length=runconfig.segment_length, 
                                       segment_trim=True, 
                                       trim5=runconfig.trim5,
                                       trim3=runconfig.trim3)
    segmented_read_length = segments[-1][1]
    logging.info("Segmented alignment will use effective read length of %d" % 
                 (segmented_read_length))
    if all(up_to_date(discordant_bam_file, fq) for fq in runconfig.fastq_files):
        logging.info("[SKIPPED] Discordant alignment results exist")
    else:
        logging.info("Aligning initially unmapped reads in single read mode")
        align(unaligned_fastq_files, runconfig.fastq_format, bowtie_index,
              discordant_bam_file, 
              bowtie_bin=runconfig.bowtie_bin,
              num_processors=runconfig.num_processors, 
              segment_length=runconfig.segment_length,
              segment_trim=True,
              trim5=runconfig.trim5, 
              trim3=runconfig.trim3, 
              multihits=runconfig.multihits,
              mismatches=runconfig.mismatches, 
              bowtie_mode=bowtie_mode)
    #
    # Merge paired-end reads step
    #
    paired_bam_file = os.path.join(runconfig.output_dir, config.DISCORDANT_PAIRED_BAM_FILE)
    if up_to_date(paired_bam_file, discordant_bam_file):
        logging.info("[SKIPPED] Read pairing results exist")
    else:
        logging.info("Pairing aligned reads")
        bamfh = pysam.Samfile(discordant_bam_file, "rb")
        paired_bamfh = pysam.Samfile(paired_bam_file, "wb", template=bamfh)
        merge_read_pairs(bamfh, paired_bamfh, 
                         runconfig.min_fragment_length,
                         runconfig.max_fragment_length,
                         library_type)
        paired_bamfh.close() 
        bamfh.close()
    #
    # Find discordant reads step
    #
    discordant_bedpe_file = os.path.join(tmp_dir, config.DISCORDANT_BEDPE_FILE)
    padding = original_read_length - segmented_read_length
    if (up_to_date(discordant_bedpe_file, paired_bam_file)):
        logging.info("[SKIPPED] Finding discordant reads")
    else:
        logging.info("Finding discordant reads")
        find_discordant_reads(pysam.Samfile(paired_bam_file, "rb"), 
                              discordant_bedpe_file, gene_feature_file,
                              max_indel_size=runconfig.max_indel_size,
                              max_isize=runconfig.max_fragment_length,
                              max_multihits=runconfig.multihits,
                              library_type=library_type,
                              padding=padding)
    #
    # Extract full sequences of the discordant reads
    #
    extended_discordant_bedpe_file = os.path.join(runconfig.output_dir, config.EXTENDED_DISCORDANT_BEDPE_FILE)
    if up_to_date(extended_discordant_bedpe_file, discordant_bedpe_file):
        logging.info("[SKIPPED] Retrieving full length sequences for realignment")
    else:
        logging.info("Retrieving full length sequences for realignment")
        extend_sequences(unaligned_fastq_files, 
                         discordant_bedpe_file,
                         extended_discordant_bedpe_file)
    #
    # Sort discordant reads
    #
    sorted_discordant_bedpe_file = os.path.join(runconfig.output_dir, config.SORTED_DISCORDANT_BEDPE_FILE)
    if (up_to_date(sorted_discordant_bedpe_file, extended_discordant_bedpe_file)):
        logging.info("[SKIPPED] Sorting discordant BEDPE file")
    else:        
        logging.info("Sorting discordant BEDPE file")
        sort_discordant_reads(extended_discordant_bedpe_file, sorted_discordant_bedpe_file)        
    #
    # Nominate chimeras step
    #
    encompassing_bedpe_file = os.path.join(tmp_dir, config.ENCOMPASSING_CHIMERA_BEDPE_FILE)        
    if (up_to_date(encompassing_bedpe_file, sorted_discordant_bedpe_file)):
        logging.info("[SKIPPED] Nominating chimeras from discordant reads")
    else:        
        logging.info("Nominating chimeras from discordant reads")
        nominate_chimeras(open(sorted_discordant_bedpe_file, "r"),
                          open(encompassing_bedpe_file, "w"),
                          gene_feature_file,                          
                          trim=config.EXON_JUNCTION_TRIM_BP)
    #
    # Filter encompassing chimeras step
    #
    filtered_encomp_bedpe_file = \
        os.path.join(runconfig.output_dir,
                     config.FILTERED_ENCOMPASSING_CHIMERA_BEDPE_FILE)
    if (up_to_date(filtered_encomp_bedpe_file, encompassing_bedpe_file)):
        logging.info("[SKIPPED] Filtering encompassing chimeras")
    else:
        logging.info("Filtering encompassing chimeras")
        # add standard deviations above the mean
        max_isize = isize_mean + runconfig.filter_isize_stdevs*isize_std
        filter_encompassing_chimeras(encompassing_bedpe_file,
                                     filtered_encomp_bedpe_file,
                                     gene_feature_file,
                                     max_multimap=runconfig.filter_max_multimaps,
                                     multimap_cov_ratio=runconfig.filter_multimap_ratio,
                                     max_isize=max_isize,
                                     strand_pval=runconfig.filter_strand_pval)
    #
    # Nominate spanning reads step
    #
    spanning_fastq_file = os.path.join(runconfig.output_dir, config.SPANNING_FASTQ_FILE)
    if (up_to_date(spanning_fastq_file, extended_discordant_bedpe_file) and 
        up_to_date(spanning_fastq_file, filtered_encomp_bedpe_file)):
        logging.info("[SKIPPED] Nominating junction spanning reads")
    else:
        logging.info("Nominating junction spanning reads")
        nominate_spanning_reads(open(extended_discordant_bedpe_file, 'r'),
                                open(filtered_encomp_bedpe_file, 'r'),
                                open(spanning_fastq_file, 'w'))    
    #
    # Extract junction sequences from chimeras file
    #        
    ref_fasta_file = os.path.join(runconfig.index_dir, config.ALIGN_INDEX + ".fa")
    junc_fasta_file = os.path.join(tmp_dir, config.JUNC_REF_FASTA_FILE)
    junc_map_file = os.path.join(tmp_dir, config.JUNC_REF_MAP_FILE)
    spanning_read_length = get_read_length(spanning_fastq_file)    
    if (up_to_date(junc_fasta_file, filtered_encomp_bedpe_file) and
        up_to_date(junc_map_file, filtered_encomp_bedpe_file)):        
        logging.info("[SKIPPED] Extracting junction read sequences")
    else:        
        logging.info("Extracting junction read sequences")
        bedpe_to_junction_fasta(filtered_encomp_bedpe_file, ref_fasta_file,                                
                                spanning_read_length, 
                                open(junc_fasta_file, "w"),
                                open(junc_map_file, "w"))
    #
    # Build a bowtie index to align and detect spanning reads
    #
    bowtie_spanning_index = os.path.join(tmp_dir, config.JUNC_BOWTIE_INDEX)
    bowtie_spanning_index_file = os.path.join(tmp_dir, config.JUNC_BOWTIE_INDEX_FILE)
    if (up_to_date(bowtie_spanning_index_file, junc_fasta_file)):
        logging.info("[SKIPPED] Building bowtie index for junction-spanning reads")
    else:        
        logging.info("Building bowtie index for junction-spanning reads")
        args = [runconfig.bowtie_build_bin, junc_fasta_file, bowtie_spanning_index]
        f = open(os.path.join(log_dir, "bowtie_build.log"), "w")
        subprocess.call(args, stdout=f, stderr=f)
        f.close()
    #
    # Align unmapped reads across putative junctions
    #
    junc_bam_file = os.path.join(runconfig.output_dir, config.JUNC_READS_BAM_FILE)
    if (up_to_date(junc_bam_file, bowtie_spanning_index_file) and
        up_to_date(junc_bam_file, spanning_fastq_file)):
        logging.info("[SKIPPED] Aligning junction spanning reads")
    else:            
        logging.info("Aligning junction spanning reads")
        retcode = align_sr_full(spanning_fastq_file, 
                                bowtie_spanning_index,
                                junc_bam_file,
                                trim5=runconfig.trim5,
                                trim3=runconfig.trim3,                                 
                                num_processors=runconfig.num_processors,
                                fastq_format=runconfig.fastq_format,
                                multihits=runconfig.multihits,
                                mismatches=runconfig.mismatches,
                                bowtie_bin=runconfig.bowtie_bin,
                                bowtie_mode=bowtie_mode)
        if retcode != 0:
            logging.error("Bowtie failed with error code %d" % (retcode))    
            sys.exit(retcode)
    #
    # Merge spanning and encompassing read information
    #
    raw_chimera_bedpe_file = os.path.join(runconfig.output_dir, config.RAW_CHIMERA_BEDPE_FILE)
    if (up_to_date(raw_chimera_bedpe_file, junc_bam_file) and
        up_to_date(raw_chimera_bedpe_file, junc_map_file)):
        logging.info("[SKIPPED] Merging spanning and encompassing read alignments")
    else:
        logging.info("Merging spanning and encompassing read alignments")
        merge_spanning_alignments(junc_bam_file, junc_map_file, 
                                  raw_chimera_bedpe_file,
                                  spanning_read_length, 
                                  anchor_min=0, 
                                  anchor_max=0,
                                  anchor_mismatches=0)
    #
    # Final filters
    #
    chimera_bedpe_file = os.path.join(runconfig.output_dir, config.CHIMERA_BEDPE_FILE)
    if (up_to_date(chimera_bedpe_file, raw_chimera_bedpe_file)):
        logging.info("[SKIPPED] Filtering chimeras")
    else:
        logging.info("Filtering chimeras")
        # add standard deviations above the mean
        max_isize = isize_mean + config.ISIZE_NUM_STDEVS*isize_std
        filter_spanning_chimeras(raw_chimera_bedpe_file, 
                                 chimera_bedpe_file,
                                 gene_feature_file)
    return JOB_SUCCESS

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    # parse run parameters in config file and command line
    runconfig = RunConfig()
    runconfig.from_args(sys.argv[1:])
    # run chimerascan
    sys.exit(run_chimerascan(runconfig))

if __name__ == '__main__':
    main()
