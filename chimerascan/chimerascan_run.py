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
from chimerascan import __version__

__author__ = "Matthew Iyer"
__copyright__ = "Copyright 2011, chimerascan project"
__credits__ = ["Matthew Iyer", "Christopher Maher"]
__license__ = "GPL"
__maintainer__ = "Matthew Iyer"
__email__ = "mkiyer@med.umich.edu"
__status__ = "beta"

import logging
import os
import subprocess
import sys
import shutil
from optparse import OptionParser, OptionGroup
import xml.etree.ElementTree as etree 

# check for python version 2.6.0 or greater
if sys.version_info < (2,6,0):
    sys.stderr.write("You need python 2.6 or later to run chimerascan\n")
    sys.exit(1)

# local imports
from chimerascan import pysam
import chimerascan.lib.config as config
from chimerascan.lib.config import JOB_SUCCESS, JOB_ERROR, MIN_SEGMENT_LENGTH
from chimerascan.lib.base import LibraryTypes, check_executable, \
    get_read_length, parse_bool, indent_xml, up_to_date
from chimerascan.lib.seq import FASTQ_QUAL_FORMATS, SANGER_FORMAT
from chimerascan.lib.fragment_size_distribution import InsertSizeDistribution

from chimerascan.pipeline.fastq_inspect_reads import inspect_reads
from chimerascan.pipeline.align_bowtie import align_pe, align_sr, trim_align_pe_sr
from chimerascan.pipeline.find_discordant_reads import find_discordant_fragments
from chimerascan.pipeline.discordant_reads_to_bedpe import discordant_reads_to_bedpe, sort_bedpe
from chimerascan.pipeline.nominate_chimeras import nominate_chimeras
from chimerascan.pipeline.chimeras_to_breakpoints import chimeras_to_breakpoints

from chimerascan.pipeline.nominate_spanning_reads import nominate_encomp_spanning_reads, nominate_unmapped_spanning_reads
from chimerascan.pipeline.merge_spanning_alignments import merge_spanning_alignments
from chimerascan.pipeline.filter_chimeras import filter_chimeras
from chimerascan.pipeline.write_output import write_output

# defaults for bowtie
DEFAULT_NUM_PROCESSORS = config.BASE_PROCESSORS
DEFAULT_BOWTIE_PATH = ""
DEFAULT_BOWTIE_ARGS = "--best --strata"
DEFAULT_DISCORD_BOWTIE_ARGS = "--best"
DEFAULT_MULTIHITS = 100
DEFAULT_MISMATCHES = 2
DEFAULT_DISCORD_MISMATCHES = 3
DEFAULT_SEGMENT_LENGTH = 25
DEFAULT_TRIM5 = 0
DEFAULT_TRIM3 = 0
DEFAULT_MIN_FRAG_LENGTH = 0
DEFAULT_MAX_FRAG_LENGTH = 1000 
DEFAULT_FASTQ_QUAL_FORMAT = SANGER_FORMAT
DEFAULT_LIBRARY_TYPE = LibraryTypes.FR_UNSTRANDED

DEFAULT_ISIZE_MEAN = 200
DEFAULT_ISIZE_STDEV = 40
DEFAULT_HOMOLOGY_MISMATCHES = config.BREAKPOINT_HOMOLOGY_MISMATCHES
DEFAULT_ANCHOR_MIN = 4
DEFAULT_ANCHOR_LENGTH = 8
DEFAULT_ANCHOR_MISMATCHES = 0
DEFAULT_FILTER_ISIZE_PERCENTILE = 99.0
DEFAULT_FILTER_UNIQUE_FRAGS = 2.0
DEFAULT_FILTER_ISOFORM_FRACTION = 0.10
NUM_POSITIONAL_ARGS = 4
DEFAULT_KEEP_TMP = True

def check_fastq_files(parser, fastq_files, segment_length, trim5, trim3):
    """
    Assures FASTQ files exist and have long enough reads to be
    trimmed by 'trim5' and 'trim3' bases on 5' and 3' ends, respectively,
    to produce reads of 'segment_length'
    """
    # check that input fastq files exist
    read_lengths = []
    for mate,fastq_file in enumerate(fastq_files):
        if not os.path.isfile(fastq_file):
            parser.error("mate '%d' fastq file '%s' is not valid" % 
                         (mate, fastq_file))
        logging.debug("Checking read length for file %s" % (fastq_file))
        read_lengths.append(get_read_length(fastq_file))
        logging.debug("Read length for file %s: %d" % 
                      (fastq_file, read_lengths[-1]))
    # check that mate read lengths are equal
    if len(set(read_lengths)) > 1:
        logging.error("read lengths mate1=%d and mate2=%d are unequal" % 
                      (read_lengths[0], read_lengths[1]))
        return False
    rlen = read_lengths[0]
    trimmed_rlen = rlen - trim5 - trim3
    # check that segment length >= MIN_SEGMENT_LENGTH 
    if segment_length < MIN_SEGMENT_LENGTH:
        parser.error("segment length (%d) too small (min is %d)" % 
                     (segment_length, MIN_SEGMENT_LENGTH))
    # check that segment length < trimmed read length
    if segment_length > trimmed_rlen:
        parser.error("segment length (%d) longer than trimmed read length (%d)" % 
                     (segment_length, trimmed_rlen))
    logging.debug("Segmented alignment will use effective read length of %d" % 
                  (segment_length))        
    return read_lengths[0]

def check_command_line_args(options, args, parser):
    # check command line arguments
    if len(args) < NUM_POSITIONAL_ARGS:
        parser.error("Incorrect number of command line arguments")
    index_dir = args[0]
    fastq_files = args[1:3]
    output_dir = args[3]
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
    # perform additional checks to ensure fastq files can be segmented
    check_fastq_files(parser, fastq_files, options.segment_length, 
                      options.trim5, options.trim3)
    # check that output dir is not a regular file
    if os.path.exists(output_dir) and (not os.path.isdir(output_dir)):
        parser.error("Output directory name '%s' exists and is not a valid directory" % 
                     (output_dir))
    if check_executable(os.path.join(options.bowtie_path, "bowtie-build")):
        logging.debug("Checking for 'bowtie-build' binary... found")
    else:
        parser.error("bowtie-build binary not found or not executable")
    # check that bowtie program exists
    if check_executable(os.path.join(options.bowtie_path, "bowtie")):
        logging.debug("Checking for 'bowtie' binary... found")
    else:
        parser.error("bowtie binary not found or not executable")
    # check that alignment index exists
    if os.path.isdir(index_dir):
        logging.debug("Checking for chimerascan index directory... found")
    else:
        parser.error("chimerascan alignment index directory '%s' not valid" % 
                     (index_dir))
    # check that alignment index file exists
    align_index_file = os.path.join(index_dir, config.BOWTIE_INDEX_FILE)
    if os.path.isfile(align_index_file):
        logging.debug("Checking for bowtie index file... found")
    else:
        parser.error("chimerascan bowtie index file '%s' invalid" % (align_index_file))
    # check for sufficient processors
    if options.num_processors < config.BASE_PROCESSORS:
        logging.warning("Please specify >=2 processes using '-p' to allow program to run efficiently")

class RunConfig(object):
    
    attrs = (("num_processors", int, DEFAULT_NUM_PROCESSORS),
             ("quals", str, DEFAULT_FASTQ_QUAL_FORMAT),
             ("keep_tmp", parse_bool, DEFAULT_KEEP_TMP),
             ("bowtie_path", str, DEFAULT_BOWTIE_PATH),
             ("bowtie_args", str, DEFAULT_BOWTIE_ARGS),
             ("discord_bowtie_args", str, DEFAULT_DISCORD_BOWTIE_ARGS),             
             ("multihits", int, DEFAULT_MULTIHITS),
             ("mismatches", int, DEFAULT_MISMATCHES),
             ("discord_mismatches", int, DEFAULT_DISCORD_MISMATCHES),
             ("segment_length", int, DEFAULT_SEGMENT_LENGTH),
             ("trim5", int, DEFAULT_TRIM5),
             ("trim3", int, DEFAULT_TRIM3),
             ("min_fragment_length", int, DEFAULT_MIN_FRAG_LENGTH),
             ("max_fragment_length", int, DEFAULT_MAX_FRAG_LENGTH),
             ("library_type", str, DEFAULT_LIBRARY_TYPE),
             ("isize_mean", int, DEFAULT_ISIZE_MEAN),
             ("isize_stdev", float, DEFAULT_ISIZE_STDEV),
             ("homology_mismatches", int, DEFAULT_HOMOLOGY_MISMATCHES),
             ("anchor_min", int, DEFAULT_ANCHOR_MIN),
             ("anchor_length", int, DEFAULT_ANCHOR_LENGTH),
             ("anchor_mismatches", int, DEFAULT_ANCHOR_MISMATCHES),
             ("filter_unique_frags", float, DEFAULT_FILTER_UNIQUE_FRAGS),
             ("filter_isize_percentile", float, DEFAULT_FILTER_ISIZE_PERCENTILE),
             ("filter_isoform_fraction", float, DEFAULT_FILTER_ISOFORM_FRACTION),
             ("filter_false_pos_file", float, ""))

    def __init__(self):
        self.output_dir = None
        self.fastq_files = None
        self.index_dir = None
        for attrname, attrtype, attrdefault in self.attrs:
            setattr(self, attrname, None)

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
        for attrname, attrtype, attrdefault in self.attrs:
            val = root.findtext(attrname, attrdefault)            
            setattr(self, attrname, attrtype(val))
    
    def to_xml(self):
        root = etree.Element("chimerascan_config")
        # output dir
        elem = etree.SubElement(root, "output_dir")
        elem.text = self.output_dir   
        # fastq files
        elem = etree.SubElement(root, "fastq_files")
        for mate,fastq_file in enumerate(self.fastq_files):
            file_elem = etree.SubElement(elem, "file", mate=str(mate))
            file_elem.text = fastq_file        
        # index
        elem = etree.SubElement(root, "index")
        elem.text = self.index_dir
        # optional arguments
        for attrname, attrtype, attrdefault in self.attrs:
            val = getattr(self, attrname)
            elem = etree.SubElement(root, attrname)
            elem.text = str(val)
        # indent for pretty printing
        indent_xml(root)        
        return etree.tostring(root)

    @staticmethod
    def get_option_parser():
        parser = OptionParser(usage="%prog [options] [--config <config_file> "
                              " | <index> <mate1.fq> <mate2.fq> <output_dir>]",
                              version="%s" % __version__)
        # standard options
        parser.add_option("--config", dest="config_file",
                          help="Path to configuration XML file") 
        parser.add_option("-v", "--verbose", dest="verbose",
                          action="store_true", default=False,
                          help="enable verbose logging output "
                          "[default=%default]")
        parser.add_option("-p", "--processors", dest="num_processors", 
                          type="int", default=DEFAULT_NUM_PROCESSORS,
                          help="Number of processor cores to allocate to "
                          "chimerascan [default=%default]")
        parser.add_option("--keep-tmp", dest="keep_tmp", action="store_true",
                          default=DEFAULT_KEEP_TMP,
                          help="DO NOT delete intermediate files after run")
        parser.add_option("--rm-tmp", dest="keep_tmp", action="store_false", 
                          help="Delete intermediate files after run")
        parser.add_option("--quals", dest="quals",
                          choices=FASTQ_QUAL_FORMATS, 
                          default=DEFAULT_FASTQ_QUAL_FORMAT, metavar="FMT",
                          help="FASTQ quality score format [default=%default]")
        parser.add_option('--library-type', dest="library_type", 
                          choices=LibraryTypes.choices(),
                          default=DEFAULT_LIBRARY_TYPE,
                          help="Library type [default=%default]")
        parser.add_option("--isize-mean", dest="isize_mean", type="int",
                          default=DEFAULT_ISIZE_MEAN, metavar="N",
                          help="Mean insert size to sample from when "
                          "insert size distribution cannot be determined "
                          "empirically [default=%default]")
        parser.add_option("--isize-stdev", dest="isize_stdev", type="float",
                          default=DEFAULT_ISIZE_STDEV, metavar="N",
                          help="Insert size standard deviation to sample "
                          "from when insert size distribution cannot be "
                          "determined empirically [default=%default]")
        # alignment options
        bowtie_group = OptionGroup(parser, "Bowtie options",
                                   "Adjust these options to change "
                                   "bowtie alignment settings")
        bowtie_group.add_option("--bowtie-path",
                                dest="bowtie_path",
                                default=DEFAULT_BOWTIE_PATH,
                                help="Path to directory containing 'bowtie' "
                                "[default: assumes bowtie is in the current "
                                "PATH]")
        bowtie_group.add_option("--trim5", type="int", dest="trim5", 
                                default=DEFAULT_TRIM5, metavar="N",
                                help="Trim N bases from 5' end of read")
        bowtie_group.add_option("--trim3", type="int", dest="trim3", 
                                default=DEFAULT_TRIM3, metavar="N",
                                help="Trim N bases from 3' end of read")
        bowtie_group.add_option("--min-fragment-length", type="int", 
                                dest="min_fragment_length", 
                                default=DEFAULT_MIN_FRAG_LENGTH,
                                help="Smallest expected fragment length "
                                "[default=%default]")
        bowtie_group.add_option("--max-fragment-length", type="int", 
                                dest="max_fragment_length", 
                                default=DEFAULT_MAX_FRAG_LENGTH,
                                help="Largest expected fragment length (reads less"
                                " than this fragment length are assumed to be"
                                " unspliced and contiguous) [default=%default]")
        bowtie_group.add_option("--multihits", type="int", dest="multihits", 
                                default=DEFAULT_MULTIHITS, metavar="N",
                                help="Ignore reads that map to more than N "
                                "locations [default=%default]")
        bowtie_group.add_option("--initial-mismatches", type="int", dest="mismatches",
                                default=DEFAULT_MISMATCHES, metavar="N",
                                help="Aligned reads must have <= N mismatches "
                                "[default=%default]")
        bowtie_group.add_option("--initial-bowtie-args", dest="bowtie_args",
                                default=DEFAULT_BOWTIE_ARGS,
                                help="Additional arguments to pass to bowtie "
                                "aligner (see Bowtie manual for options) "
                                "[default='%default']")
        bowtie_group.add_option("--discord-bowtie-args", dest="discord_bowtie_args",
                                default=DEFAULT_DISCORD_BOWTIE_ARGS,
                                help="Additional arguments to pass to bowtie "
                                "aligner for discordant alignment phase (see "
                                "Bowtie manual for options) [default='%default']")
        bowtie_group.add_option("--discord-mismatches", type="int", dest="discord_mismatches",
                                default=DEFAULT_DISCORD_MISMATCHES, metavar="N",
                                help="Discordant aligned reads must have <= N mismatches "
                                "[default=%default]")
        bowtie_group.add_option("--segment-length", type="int", 
                                dest="segment_length", 
                                default=DEFAULT_SEGMENT_LENGTH,
                                metavar="N",
                                help="Size of read segments during discordant " 
                                "alignment phase [default=%default]")
        parser.add_option_group(bowtie_group)
        # filtering options
        filter_group = OptionGroup(parser, "Filtering options",
                                   "Adjust these options to change "
                                   "filtering behavior")
        filter_group.add_option("--homology-mismatches", type="int", 
                                dest="homology_mismatches",
                                default=DEFAULT_HOMOLOGY_MISMATCHES,
                                help="Number of mismatches to tolerate at "
                                "breakpoints when computing homology "
                                "[default=%default]", metavar="N")
        filter_group.add_option("--anchor-min", type="int", dest="anchor_min", 
                                default=DEFAULT_ANCHOR_MIN, metavar="N",
                                help="Minimum breakpoint overlap (bp) required "
                                "to call spanning reads [default=%default]")
        filter_group.add_option("--anchor-length", type="int", dest="anchor_length", 
                                default=DEFAULT_ANCHOR_LENGTH, metavar="N",
                                help="Size (bp) of anchor region where mismatch "
                                "checks are enforced [default=%default]")
        filter_group.add_option("--anchor-mismatches", type="int", dest="anchor_mismatches", 
                                default=DEFAULT_ANCHOR_MISMATCHES,
                                metavar="N",
                                help="Number of mismatches allowed within anchor "
                                "region [default=%default]")
        filter_group.add_option("--filter-unique-frags", type="float",
                                default=DEFAULT_FILTER_UNIQUE_FRAGS,
                                dest="filter_unique_frags", metavar="N",
                                help="Filter chimeras with less than N unique "
                                "aligned fragments (multimapping fragments are "
                                "assigned fractional weights) [default=%default]")
        filter_group.add_option("--filter-isize-percentile", type="float",
                                default=DEFAULT_FILTER_ISIZE_PERCENTILE,
                                dest="filter_isize_percentile", metavar="N",
                                help="Filter chimeras when putative insert "
                                "size is larger than the Nth percentile "
                                "of the distribution [default=%default]")
        filter_group.add_option("--filter-isoform-fraction", type="float", 
                                default=DEFAULT_FILTER_ISOFORM_FRACTION, metavar="X",
                                help="Filter chimeras with expression ratio "
                                " less than X (0.0-1.0) relative to the wild-type "
                                "5' transcript level [default=%default]")
        filter_group.add_option("--filter-false-pos", default="",
                                dest="filter_false_pos_file",
                                help="File containing known false positive "
                                "chimeric transcript pairs to filter out")
        parser.add_option_group(filter_group)
        return parser        

    def from_args(self, args, parser=None):
        if parser is None:
            parser = self.get_option_parser()
        options, args = parser.parse_args(args=args)
        # parse config file options/args
        if options.config_file is not None:
            self.from_xml(options.config_file)
        # reset logging to verbose
        if options.verbose:
            logging.getLogger().setLevel(logging.DEBUG)
        # check command line arguments
        if self.index_dir is None:
            self.index_dir = args[0]            
        if self.fastq_files is None:
            if len(args) < 2:
                parser.error("fastq files not specified in config file or command line")        
            self.fastq_files = args[1:3]
        if self.output_dir is None:
            if len(args) < 3:
                parser.error("output dir not specified in config file or command line")   
            self.output_dir = args[3]
        # optional arguments
        # set rest of options, overriding if attribute is undefined
        # or set to something other than the default 
        for attrname, attrtype, attrdefault in self.attrs:
            if ((getattr(self, attrname) is None) or
                (getattr(options, attrname) != attrdefault)):
                setattr(self, attrname, getattr(options, attrname))        

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
        if check_executable(os.path.join(self.bowtie_path, "bowtie-build")):
            logging.debug("Checking for 'bowtie-build' binary... found")
        else:
            logging.error("bowtie-build binary not found or not executable")
            config_passed = False
        # check that bowtie program exists
        if check_executable(os.path.join(self.bowtie_path, "bowtie")):
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

def run_chimerascan(runconfig):
    """
    main function for running the chimerascan pipeline
    """
    # print a welcome message
    title_string = "Running chimerascan version %s" % (__version__)
    logging.info(title_string)
    logging.info("-" * len(title_string))
    # validate run configuration
    config_passed = runconfig.check_config()
    if not config_passed:
        logging.error("Invalid run configuration, aborting.")
        return JOB_ERROR
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
    # write the run config to a file
    xmlstring = runconfig.to_xml()
    runconfig_xml_file = os.path.join(runconfig.output_dir, config.RUNCONFIG_XML_FILE)
    logging.info("Writing run configuration to XML file: %s" % (runconfig_xml_file))
    fh = open(runconfig_xml_file, "w")
    print >>fh, xmlstring
    fh.close()
    # gather and parse run parameters
    bowtie_index = os.path.join(runconfig.index_dir, config.ALIGN_INDEX)
    bowtie_bin = os.path.join(runconfig.bowtie_path, "bowtie")
    bowtie_build_bin = os.path.join(runconfig.bowtie_path, "bowtie-build")
    original_read_length = get_read_length(runconfig.fastq_files[0])
    # minimum fragment length cannot be smaller than the trimmed read length
    trimmed_read_length = original_read_length - runconfig.trim5 - runconfig.trim3
    min_fragment_length = max(runconfig.min_fragment_length, 
                              trimmed_read_length)
    # 
    # Process and inspect the FASTQ files, performing several alterations 
    # to the reads:
    #
    # 1) rename them from long string to numbers to save space throughout
    # the pipeline
    # 2) ensure the "/1" and "/2" suffixes exist to denote reads in a pair
    # 3) convert quality scores to sanger format
    # 4) #TODO# trim poor quality reads 
    # (requires support for variable length reads)
    # 
    converted_fastq_files = [os.path.join(tmp_dir, fq) 
                             for fq in config.CONVERTED_FASTQ_FILES]
    msg = "Processing FASTQ files"
    if (all(up_to_date(cfq, fq) for cfq,fq in 
            zip(converted_fastq_files, runconfig.fastq_files))):
        logging.info("[SKIPPED] %s" % (msg))
    else:
        logging.info(msg)
        converted_fastq_prefix = \
            os.path.join(tmp_dir, config.CONVERTED_FASTQ_PREFIX)
        try:
            retcode = inspect_reads(runconfig.fastq_files, 
                                    converted_fastq_prefix,
                                    quals=runconfig.quals)
            if retcode != config.JOB_SUCCESS:
                logging.error("%s step failed" % (msg))
                return JOB_ERROR
        except Exception as e:
            logging.info("Cleaning up after error %s" % (str(e)))
            for fq in converted_fastq_files:
                if os.path.isfile(fq):
                    os.remove(fq)
    #
    # Initial Bowtie alignment step
    #
    # align in paired-end mode, trying to resolve as many reads as possible
    # this effectively rules out the vast majority of reads as candidate
    # fusions
    #
    unaligned_fastq_param = os.path.join(tmp_dir, config.UNALIGNED_FASTQ_PARAM)
    unaligned_fastq_files = [os.path.join(tmp_dir, fq) 
                             for fq in config.UNALIGNED_FASTQ_FILES]
    maxmultimap_fastq_param = os.path.join(tmp_dir, config.MAXMULTIMAP_FASTQ_PARAM)
    aligned_bam_file = os.path.join(runconfig.output_dir, config.ALIGNED_READS_BAM_FILE)
    aligned_log_file = os.path.join(log_dir, "bowtie_alignment.log")    
    msg = "Aligning reads with bowtie"
    if (all(up_to_date(aligned_bam_file, fq) for fq in converted_fastq_files) and
        all(up_to_date(a,b) for a,b in zip(unaligned_fastq_files, converted_fastq_files))):
        logging.info("[SKIPPED] %s" % (msg))
    else:    
        logging.info(msg)
        retcode = align_pe(converted_fastq_files, 
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
                           quals=SANGER_FORMAT,
                           multihits=runconfig.multihits,
                           mismatches=runconfig.mismatches,
                           bowtie_bin=bowtie_bin,
                           bowtie_args=runconfig.bowtie_args,
                           log_file=aligned_log_file,
                           keep_unmapped=False)
        if retcode != 0:
            logging.error("Bowtie failed with error code %d" % (retcode))
            return JOB_ERROR
    #
    # Sort aligned reads by position
    #
    msg = "Sorting aligned reads"
    sorted_aligned_bam_file = os.path.join(runconfig.output_dir, 
                                           config.SORTED_ALIGNED_READS_BAM_FILE)
    if (up_to_date(sorted_aligned_bam_file, aligned_bam_file)):
        logging.info("[SKIPPED] %s" % (msg))
    else:
        logging.info(msg)
        sorted_aligned_bam_prefix = os.path.splitext(sorted_aligned_bam_file)[0]
        pysam.sort("-m", str(int(1e9)), aligned_bam_file, sorted_aligned_bam_prefix)
    #
    # Index BAM file
    #
    msg = "Indexing BAM file"
    sorted_aligned_bam_index_file = sorted_aligned_bam_file + ".bai"
    if (up_to_date(sorted_aligned_bam_index_file, sorted_aligned_bam_file)):
        logging.info("[SKIPPED] %s" % (msg))
    else:
        logging.info(msg)
        pysam.index(sorted_aligned_bam_file)
    #
    # Get insert size distribution
    #
    isize_dist_file = os.path.join(runconfig.output_dir, 
                                   config.ISIZE_DIST_FILE)
    if up_to_date(isize_dist_file, aligned_bam_file):
        logging.info("[SKIPPED] Profiling insert size distribution")
        isize_dist = InsertSizeDistribution.from_file(open(isize_dist_file, "r"))
    else:
        logging.info("Profiling insert size distribution")
        max_isize_samples = config.ISIZE_MAX_SAMPLES
        bamfh = pysam.Samfile(aligned_bam_file, "rb")
        isize_dist = InsertSizeDistribution.from_bam(bamfh, min_isize=min_fragment_length, 
                                                     max_isize=runconfig.max_fragment_length, 
                                                     max_samples=max_isize_samples)
        bamfh.close()
        # if not enough samples, use a normal distribution instead
        # of the empirical distribution
        if isize_dist.n < config.ISIZE_MIN_SAMPLES:
            logging.warning("Not enough fragments to sample insert size "
                            "distribution empirically.  Using mean=%d "
                            "stdev=%f instead" % 
                            (runconfig.isize_mean, 
                             runconfig.isize_stdev))
            isize_dist = InsertSizeDistribution.from_random(runconfig.isize_mean, 
                                                            runconfig.isize_stdev, 
                                                            min_isize=runconfig.min_fragment_length,
                                                            max_isize=runconfig.max_fragment_length,
                                                            samples=max_isize_samples)
        isize_dist.to_file(open(isize_dist_file, "w"))
    logging.info("Insert size samples=%d mean=%f std=%f median=%d mode=%d" % 
                 (isize_dist.n, isize_dist.mean(), isize_dist.std(), 
                  isize_dist.isize_at_percentile(50.0), isize_dist.mode()))
    #
    # Realignment step
    #
    # trim and realign all the initially unaligned reads, in order to
    # increase sensitivity to detect reads spanning fusion junctions
    #
    realigned_bam_file = os.path.join(tmp_dir, config.REALIGNED_BAM_FILE)
    realigned_log_file = os.path.join(log_dir, "bowtie_trimmed_realignment.log")
    msg = "Trimming and re-aligning initially unmapped reads"
    if all(up_to_date(realigned_bam_file, fq) 
           for fq in unaligned_fastq_files):
        logging.info("[SKIPPED] %s" % (msg))
    else:
        logging.info(msg)
        trim_align_pe_sr(unaligned_fastq_files,
                         bowtie_index,
                         realigned_bam_file,
                         unaligned_fastq_param=None,
                         maxmultimap_fastq_param=None,
                         trim5=runconfig.trim5,
                         library_type=runconfig.library_type,
                         num_processors=runconfig.num_processors,
                         quals=SANGER_FORMAT,
                         multihits=runconfig.multihits, 
                         mismatches=runconfig.discord_mismatches,
                         bowtie_bin=bowtie_bin,
                         bowtie_args=runconfig.discord_bowtie_args,
                         log_file=realigned_log_file,
                         segment_length=runconfig.segment_length,
                         keep_unmapped=True)
    #
    # Find discordant reads
    #
    # iterate through realigned reads and divide them into groups of
    # concordant, discordant within a gene (isoforms), discordant
    # between different genes, and discordant in the genome
    #
    gene_paired_bam_file = os.path.join(tmp_dir, config.GENE_PAIRED_BAM_FILE)
    genome_paired_bam_file = os.path.join(tmp_dir, config.GENOME_PAIRED_BAM_FILE)
    realigned_unmapped_bam_file = os.path.join(tmp_dir, config.REALIGNED_UNMAPPED_BAM_FILE)
    realigned_complex_bam_file = os.path.join(tmp_dir, config.REALIGNED_COMPLEX_BAM_FILE)
    msg = "Classifying concordant and discordant read pairs"
    if (up_to_date(gene_paired_bam_file, realigned_bam_file) and
        up_to_date(genome_paired_bam_file, realigned_bam_file) and
        up_to_date(realigned_unmapped_bam_file, realigned_bam_file) and
        up_to_date(realigned_complex_bam_file, realigned_bam_file)):        
        logging.info("[SKIPPED] %s" % (msg))
    else:
        logging.info(msg)
        find_discordant_fragments(realigned_bam_file, 
                                  gene_paired_bam_file=gene_paired_bam_file, 
                                  genome_paired_bam_file=genome_paired_bam_file, 
                                  unmapped_bam_file=realigned_unmapped_bam_file,
                                  complex_bam_file=realigned_complex_bam_file,
                                  index_dir=runconfig.index_dir,
                                  max_isize=runconfig.max_fragment_length,
                                  library_type=runconfig.library_type) 

    #
    # Write a BEDPE file from discordant reads
    #
    discordant_bedpe_file = os.path.join(tmp_dir, config.DISCORDANT_BEDPE_FILE)
    msg = "Converting discordant reads to BEDPE format"
    if (up_to_date(discordant_bedpe_file, gene_paired_bam_file)):
        logging.info("[SKIPPED] %s" % (msg))
    else:
        logging.info(msg)
        discordant_reads_to_bedpe(runconfig.index_dir,
                                  gene_paired_bam_file,
                                  discordant_bedpe_file)
    #
    # Sort BEDPE file
    #
    sorted_discordant_bedpe_file = os.path.join(tmp_dir, config.SORTED_DISCORDANT_BEDPE_FILE)
    msg = "Sorting BEDPE file"
    if (up_to_date(sorted_discordant_bedpe_file, discordant_bedpe_file)):
        logging.info("[SKIPPED] %s" % (msg))
    else:
        logging.info(msg)
        # sort BEDPE file by paired chromosome/position
        sort_bedpe(input_file=discordant_bedpe_file,
                   output_file=sorted_discordant_bedpe_file,
                   tmp_dir=tmp_dir)
    #
    # Nominate chimeric genes
    #
    # iterate through pairs of reads and nominate chimeras
    #
    encompassing_chimera_file = os.path.join(tmp_dir, config.ENCOMPASSING_CHIMERA_FILE)
    msg = "Nominating chimeras from discordant reads"
    if (up_to_date(encompassing_chimera_file, sorted_discordant_bedpe_file)):
        logging.info("[SKIPPED] %s" % (msg))
    else:        
        logging.info(msg)
        retcode = nominate_chimeras(runconfig.index_dir,
                                    isize_dist_file=isize_dist_file,
                                    input_file=sorted_discordant_bedpe_file,
                                    output_file=encompassing_chimera_file,
                                    trim_bp=config.EXON_JUNCTION_TRIM_BP,
                                    min_frags=runconfig.filter_unique_frags,
                                    max_read_length=trimmed_read_length,
                                    homology_mismatches=runconfig.homology_mismatches)
        if retcode != config.JOB_SUCCESS:
            logging.error("Failed with error code %d" % (retcode))
            return config.JOB_ERROR
    #
    # Make a breakpoint FASTA file and map of breakpoints to chimeras
    # for use in spanning read alignment
    #
    breakpoint_chimera_file = os.path.join(tmp_dir, config.BREAKPOINT_CHIMERA_FILE)
    breakpoint_map_file = os.path.join(tmp_dir, config.BREAKPOINT_MAP_FILE)
    breakpoint_fasta_file = os.path.join(tmp_dir, config.BREAKPOINT_FASTA_FILE)
    msg = "Extracting breakpoint sequences from chimeras"
    if (up_to_date(breakpoint_chimera_file, encompassing_chimera_file) and
        up_to_date(breakpoint_map_file, encompassing_chimera_file) and
        up_to_date(breakpoint_fasta_file, encompassing_chimera_file)):
        logging.info("[SKIPPED] %s" % (msg))
    else:
        logging.info(msg)
        chimeras_to_breakpoints(input_file=encompassing_chimera_file,
                                breakpoint_sorted_chimera_file=breakpoint_chimera_file, 
                                breakpoint_map_file=breakpoint_map_file, 
                                breakpoint_fasta_file=breakpoint_fasta_file,
                                tmp_dir=tmp_dir)
    #
    # Build a bowtie index to align and detect spanning reads
    #
    breakpoint_bowtie_index = os.path.join(tmp_dir, config.BREAKPOINT_BOWTIE_INDEX)
    breakpoint_bowtie_index_file = os.path.join(tmp_dir, config.BREAKPOINT_BOWTIE_INDEX_FILE)
    msg = "Building bowtie index of breakpoint spanning sequences"
    if (up_to_date(breakpoint_bowtie_index_file, breakpoint_fasta_file)):
        logging.info("[SKIPPED] %s" % (msg))
    else:
        logging.info(msg)
        args = [bowtie_build_bin, breakpoint_fasta_file, 
                breakpoint_bowtie_index]
        f = open(os.path.join(log_dir, "breakpoint_bowtie_index.log"), "w")
        retcode = subprocess.call(args, stdout=f, stderr=f)
        f.close()
        if retcode != config.JOB_SUCCESS:
            logging.error("'bowtie-build' failed to create breakpoint index")
            if os.path.exists(breakpoint_bowtie_index_file):
                os.remove(breakpoint_bowtie_index_file)
            return config.JOB_ERROR
    #
    # Extract reads that were encompassing when trimmed but may span
    # breakpoints when extended.
    #
    encomp_spanning_fastq_file = os.path.join(tmp_dir, config.ENCOMP_SPANNING_FASTQ_FILE)
    msg = "Extracting encompassing reads that may extend past breakpoints"
    if (up_to_date(encomp_spanning_fastq_file, encompassing_chimera_file)):
        logging.info("[SKIPPED] %s" % (msg))
    else:
        logging.info(msg)    
        retcode = nominate_encomp_spanning_reads(encompassing_chimera_file, 
                                                 encomp_spanning_fastq_file)
        if retcode != config.JOB_SUCCESS:
            logging.error("Failed to extract breakpoint-spanning reads")
            if os.path.exists(encomp_spanning_fastq_file):
                os.remove(encomp_spanning_fastq_file)
            return config.JOB_ERROR
    #
    # Extract unmapped reads to align against breakpoint sequences
    #
    msg = "Extracting unmapped reads that may span breakpoints"
    singlemap_spanning_fastq_file = os.path.join(tmp_dir, config.SINGLEMAP_SPANNING_FASTQ_FILE)
    unaligned_spanning_fastq_file = os.path.join(tmp_dir, config.UNALIGNED_SPANNING_FASTQ_FILE)
    if (up_to_date(singlemap_spanning_fastq_file, encompassing_chimera_file, nzsize=False) and
        up_to_date(unaligned_spanning_fastq_file, encompassing_chimera_file, nzsize=False)):
        logging.info("[SKIPPED] %s" % (msg))
    else:
        logging.info(msg)
        retcode = nominate_unmapped_spanning_reads(encompassing_chimera_file, 
                                                   realigned_unmapped_bam_file,
                                                   singlemap_spanning_fastq_file,
                                                   unaligned_spanning_fastq_file,
                                                   library_type=runconfig.library_type,
                                                   tmp_dir=tmp_dir)
        if retcode != config.JOB_SUCCESS:
            logging.error("Failed to extract unmapped reads")
            if os.path.exists(unaligned_spanning_fastq_file):
                os.remove(unaligned_spanning_fastq_file)
            return config.JOB_ERROR
    #
    # Re-align encompassing reads that overlap breakpoint junctions
    # 
    encomp_spanning_bam_file = os.path.join(tmp_dir, config.ENCOMP_SPANNING_BAM_FILE)
    encomp_spanning_log_file = os.path.join(log_dir, "bowtie_encomp_spanning.log")
    msg = "Realigning encompassing reads to breakpoints"
    if (up_to_date(encomp_spanning_bam_file, breakpoint_bowtie_index_file) and
        up_to_date(encomp_spanning_bam_file, encomp_spanning_fastq_file)):
        logging.info("[SKIPPED] %s" % (msg))
    else:            
        logging.info(msg)
        retcode= align_sr(encomp_spanning_fastq_file, 
                          bowtie_index=breakpoint_bowtie_index,
                          output_bam_file=encomp_spanning_bam_file, 
                          unaligned_fastq_param=None,
                          maxmultimap_fastq_param=None,
                          trim5=0,
                          trim3=0,
                          library_type=runconfig.library_type,
                          num_processors=runconfig.num_processors, 
                          quals=SANGER_FORMAT,
                          multihits=runconfig.multihits,
                          mismatches=runconfig.mismatches, 
                          bowtie_bin=bowtie_bin,
                          bowtie_args=runconfig.bowtie_args,
                          log_file=encomp_spanning_log_file,
                          keep_unmapped=False)
        if retcode != config.JOB_SUCCESS:
            logging.error("Bowtie failed with error code %d" % (retcode))    
            return config.JOB_ERROR
    #
    # Sort encomp/spanning reads by reference/position
    #
    msg = "Sorting/indexing encompassing/spanning alignments"
    sorted_encomp_spanning_bam_file = os.path.join(tmp_dir, config.SORTED_ENCOMP_SPANNING_BAM_FILE)
    if (up_to_date(sorted_encomp_spanning_bam_file, encomp_spanning_bam_file)):
        logging.info("[SKIPPED] %s" % (msg))
    else:
        logging.info(msg)
        sorted_encomp_spanning_bam_prefix = os.path.splitext(sorted_encomp_spanning_bam_file)[0]
        pysam.sort("-m", str(int(1e9)), encomp_spanning_bam_file, sorted_encomp_spanning_bam_prefix)
        pysam.index(sorted_encomp_spanning_bam_file)
    #
    # Align single-mapping reads that may overlap breakpoint junctions
    #
    singlemap_spanning_bam_file = os.path.join(tmp_dir, config.SINGLEMAP_SPANNING_BAM_FILE)
    singlemap_spanning_log_file = os.path.join(log_dir, "bowtie_singlemap_spanning.log")
    msg = "Realigning single-mapping reads to breakpoints"
    if (up_to_date(singlemap_spanning_bam_file, breakpoint_bowtie_index_file) and
        up_to_date(singlemap_spanning_bam_file, singlemap_spanning_fastq_file)):
        logging.info("[SKIPPED] %s" % (msg))
    else:            
        logging.info(msg)
        retcode= align_sr(singlemap_spanning_fastq_file, 
                          bowtie_index=breakpoint_bowtie_index,
                          output_bam_file=singlemap_spanning_bam_file, 
                          unaligned_fastq_param=None,
                          maxmultimap_fastq_param=None,
                          trim5=0,
                          trim3=0,
                          library_type=runconfig.library_type,
                          num_processors=runconfig.num_processors, 
                          quals=SANGER_FORMAT,
                          multihits=runconfig.multihits,
                          mismatches=runconfig.mismatches, 
                          bowtie_bin=bowtie_bin,
                          bowtie_args=runconfig.bowtie_args,
                          log_file=singlemap_spanning_log_file,
                          keep_unmapped=False)
        if retcode != config.JOB_SUCCESS:
            logging.error("Bowtie failed with error code %d" % (retcode))    
            return config.JOB_ERROR
    #
    # Sort encomp/spanning reads by reference/position
    #
    msg = "Sorting/indexing single-mapping/spanning alignments"
    sorted_singlemap_spanning_bam_file = os.path.join(tmp_dir, config.SORTED_SINGLEMAP_SPANNING_BAM_FILE)
    if (up_to_date(sorted_singlemap_spanning_bam_file, singlemap_spanning_bam_file)):
        logging.info("[SKIPPED] %s" % (msg))
    else:
        logging.info(msg)
        sorted_singlemap_spanning_bam_prefix = os.path.splitext(sorted_singlemap_spanning_bam_file)[0]
        pysam.sort("-m", str(int(1e9)), singlemap_spanning_bam_file, sorted_singlemap_spanning_bam_prefix)
        pysam.index(sorted_singlemap_spanning_bam_file)
    #
    # Align unmapped mapping reads that may overlap breakpoint junctions
    #
    unaligned_spanning_bam_file = os.path.join(tmp_dir, config.UNALIGNED_SPANNING_BAM_FILE)
    unaligned_spanning_log_file = os.path.join(log_dir, "bowtie_unaligned_spanning.log")
    msg = "Realigning unmapped reads to breakpoints"
    if (up_to_date(unaligned_spanning_bam_file, breakpoint_bowtie_index_file) and
        up_to_date(unaligned_spanning_bam_file, unaligned_spanning_fastq_file)):
        logging.info("[SKIPPED] %s" % (msg))
    else:            
        logging.info(msg)
        retcode= align_sr(unaligned_spanning_fastq_file, 
                          bowtie_index=breakpoint_bowtie_index,
                          output_bam_file=unaligned_spanning_bam_file, 
                          unaligned_fastq_param=None,
                          maxmultimap_fastq_param=None,
                          trim5=0,
                          trim3=0,
                          library_type=runconfig.library_type,
                          num_processors=runconfig.num_processors, 
                          quals=SANGER_FORMAT,
                          multihits=runconfig.multihits,
                          mismatches=runconfig.mismatches, 
                          bowtie_bin=bowtie_bin,
                          bowtie_args=runconfig.bowtie_args,
                          log_file=unaligned_spanning_log_file,
                          keep_unmapped=False)
        if retcode != config.JOB_SUCCESS:
            logging.error("Bowtie failed with error code %d" % (retcode))    
            return config.JOB_ERROR
    #
    # Sort encomp/spanning reads by reference/position
    #
    msg = "Sorting/indexing unaligned/spanning alignments"
    sorted_unaligned_spanning_bam_file = os.path.join(tmp_dir, config.SORTED_UNALIGNED_SPANNING_BAM_FILE)
    if (up_to_date(sorted_unaligned_spanning_bam_file, unaligned_spanning_bam_file)):
        logging.info("[SKIPPED] %s" % (msg))
    else:
        logging.info(msg)
        sorted_unaligned_spanning_bam_prefix = os.path.splitext(sorted_unaligned_spanning_bam_file)[0]
        pysam.sort("-m", str(int(1e9)), unaligned_spanning_bam_file, sorted_unaligned_spanning_bam_prefix)
        pysam.index(sorted_unaligned_spanning_bam_file)
    #
    # Merge spanning read alignment information
    #
    spanning_chimera_file = os.path.join(tmp_dir, config.SPANNING_CHIMERA_FILE)
    msg = "Merging spanning read information"
    if (up_to_date(spanning_chimera_file, breakpoint_chimera_file) and
        up_to_date(spanning_chimera_file, unaligned_spanning_bam_file) and
        up_to_date(spanning_chimera_file, encomp_spanning_bam_file) and
        up_to_date(spanning_chimera_file, realigned_unmapped_bam_file)):
        logging.info("[SKIPPED] %s" % (msg))
    else:
        logging.info(msg)        
        merge_spanning_alignments(breakpoint_chimera_file=breakpoint_chimera_file,
                                  encomp_bam_file=sorted_encomp_spanning_bam_file,
                                  singlemap_bam_file=sorted_singlemap_spanning_bam_file,
                                  unaligned_bam_file=sorted_unaligned_spanning_bam_file,
                                  output_chimera_file=spanning_chimera_file,
                                  anchor_min=runconfig.anchor_min,
                                  anchor_length=runconfig.anchor_length,
                                  anchor_mismatches=runconfig.anchor_mismatches,
                                  library_type=runconfig.library_type,
                                  tmp_dir=tmp_dir)
    #
    # Resolve reads mapping to multiple chimeras
    # TODO: work in progress
    #
    pass
    #
    # Filter chimeras
    # 
    filtered_chimera_file = os.path.join(tmp_dir, config.FILTERED_CHIMERA_FILE)
    msg = "Filtering chimeras"
    if up_to_date(filtered_chimera_file, spanning_chimera_file):
        logging.info("[SKIPPED] %s" % (msg))
    else:
        logging.info(msg)
        # get insert size at prob    
        max_isize = isize_dist.isize_at_percentile(runconfig.filter_isize_percentile)
        median_isize = isize_dist.isize_at_percentile(50.0)
        filter_chimeras(input_file=spanning_chimera_file, 
                        output_file=filtered_chimera_file,
                        index_dir=runconfig.index_dir,
                        bam_file=sorted_aligned_bam_file,
                        unique_frags=runconfig.filter_unique_frags,
                        median_isize=median_isize,
                        isoform_fraction=runconfig.filter_isoform_fraction,
                        false_pos_file=runconfig.filter_false_pos_file)
    #
    # Write user-friendly output file
    # 
    chimera_output_file = os.path.join(runconfig.output_dir, config.CHIMERA_OUTPUT_FILE)
    msg = "Writing chimeras to file %s" % (chimera_output_file)
    if up_to_date(chimera_output_file, filtered_chimera_file):
        logging.info("[SKIPPED] %s" % (msg))
    else:
        logging.info(msg)
        write_output(filtered_chimera_file, chimera_output_file,
                     index_dir=runconfig.index_dir)
    #
    # Cleanup
    # 
    if not runconfig.keep_tmp:
        logging.info("Cleaning up temporary files")
        shutil.rmtree(tmp_dir)
    #
    # Done
    # 
    logging.info("Finished run.")
    return JOB_SUCCESS


def main():
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    # parse run parameters in config file and command line
    runconfig = RunConfig()
    runconfig.from_args(sys.argv[1:])
    # run chimerascan
    sys.exit(run_chimerascan(runconfig))

if __name__ == '__main__':
    main()
