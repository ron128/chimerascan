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
import argparse
import xml.etree.ElementTree as etree 

# check for python version 2.7.0 or greater
if sys.version_info < (2,7,1):
    sys.stderr.write("You need python 2.7.1 or later to run chimerascan\n")
    sys.exit(1)

# local imports
from chimerascan import pysam
import chimerascan.lib.config as config
from chimerascan.lib.base import LibraryTypes, check_executable, \
    parse_bool, indent_xml, up_to_date
from chimerascan.lib.seq import FASTQ_QUAL_FORMATS, SANGER_FORMAT, detect_read_length
from chimerascan.lib.fragment_size_distribution import InsertSizeDistribution

from chimerascan.pipeline.process_input_reads import process_input_reads
from chimerascan.pipeline.align_bowtie2 import bowtie2_align_pe, bowtie2_align_pe_sr, bowtie2_align_sr
from chimerascan.pipeline.find_discordant_reads import find_discordant_fragments
from chimerascan.pipeline.discordant_reads_to_bedpe import discordant_reads_to_bedpe, sort_bedpe
from chimerascan.pipeline.nominate_chimeras import nominate_chimeras
from chimerascan.pipeline.chimeras_to_breakpoints import chimeras_to_breakpoints
from chimerascan.pipeline.nominate_spanning_reads import nominate_encomp_spanning_reads, extract_single_mapped_reads, nominate_single_mapped_spanning_reads
from chimerascan.pipeline.merge_spanning_alignments import merge_spanning_alignments
from chimerascan.pipeline.resolve_discordant_reads import resolve_discordant_reads
from chimerascan.pipeline.filter_chimeras import filter_chimeras, filter_highest_coverage_isoforms, filter_encompassing_chimeras
from chimerascan.pipeline.filter_homologous_genes import filter_homologous_genes
from chimerascan.pipeline.write_output import write_output

# global default parameters
DEFAULT_NUM_PROCESSORS = config.BASE_PROCESSORS
DEFAULT_KEEP_TMP = True

# default sequencing data parameters
DEFAULT_FASTQ_QUAL_FORMAT = SANGER_FORMAT
DEFAULT_LIBRARY_TYPE = LibraryTypes.FR_UNSTRANDED
DEFAULT_ISIZE_MEAN = 200
DEFAULT_ISIZE_STDEV = 40
DEFAULT_FRAG_SIZE_SENSITIVITY = 1.0

# bowtie2 default arguments
DEFAULT_BOWTIE2_ARGS = ["--end-to-end",
                        "--very-sensitive",
                        "--reorder"]

# defaults for bowtie
DEFAULT_MIN_FRAG_LENGTH = 0
DEFAULT_MAX_FRAG_LENGTH = 1000 
DEFAULT_TRIM5 = 0
DEFAULT_TRIM3 = 0
DEFAULT_SEGMENT_LENGTH = None
DEFAULT_MAX_MULTIHITS = 1
DEFAULT_HOMOLOGY_MISMATCHES = config.BREAKPOINT_HOMOLOGY_MISMATCHES
DEFAULT_ANCHOR_MIN = 4
DEFAULT_ANCHOR_LENGTH = 8
DEFAULT_ANCHOR_MISMATCHES = 0
DEFAULT_FILTER_ISIZE_PROB = 0.01
DEFAULT_FILTER_UNIQUE_FRAGS = 2.0
DEFAULT_FILTER_ISOFORM_FRACTION = 0.01

class RunConfig(object):
    
    attrs = (("num_processors", int, DEFAULT_NUM_PROCESSORS),
             ("keep_tmp", parse_bool, DEFAULT_KEEP_TMP),
             ("quals", str, DEFAULT_FASTQ_QUAL_FORMAT),
             ("library_type", str, DEFAULT_LIBRARY_TYPE),
             ("isize_mean", int, DEFAULT_ISIZE_MEAN),
             ("isize_stdev", float, DEFAULT_ISIZE_STDEV),
             ("min_fragment_length", int, DEFAULT_MIN_FRAG_LENGTH),
             ("max_fragment_length", int, DEFAULT_MAX_FRAG_LENGTH),
             ("trim5", int, DEFAULT_TRIM5),
             ("trim3", int, DEFAULT_TRIM3),
             ("segment_length", int, DEFAULT_SEGMENT_LENGTH),
             ("max_multihits", int, DEFAULT_MAX_MULTIHITS),
             ("homology_mismatches", int, DEFAULT_HOMOLOGY_MISMATCHES),
             ("anchor_min", int, DEFAULT_ANCHOR_MIN),
             ("anchor_length", int, DEFAULT_ANCHOR_LENGTH),
             ("anchor_mismatches", int, DEFAULT_ANCHOR_MISMATCHES),
             ("filter_unique_frags", float, DEFAULT_FILTER_UNIQUE_FRAGS),
             ("filter_isize_prob", float, DEFAULT_FILTER_ISIZE_PROB),
             ("filter_isoform_fraction", float, DEFAULT_FILTER_ISOFORM_FRACTION),
             ("filter_false_pos_file", float, ""))

    def __init__(self):
        self.output_dir = None
        self.fastq_files = [None,None]
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
        root = etree.Element("chimerascan")
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
    def get_argument_parser():
        parser = argparse.ArgumentParser(usage="%(prog)s [options] <index> "
                                         "<mate1.fq> <mate2.fq> <output_dir>")
        # required options
        parser.add_argument("index_dir", default=None,
                            help="Location of chimerascan index directory")
        parser.add_argument("read1", default=None,
                            help="Path to read1 FASTQ file")
        parser.add_argument("read2", default=None,
                            help="Path to read2 FASTQ file")
        parser.add_argument("output_dir", default=None,
                            help="Location of output files")
        # standard options
        parser.add_argument('--version', action='version', 
                            version='%s' % __version__)
        parser.add_argument("--config-file", dest="config_file", 
                            help="Load parameters from a XML file "
                            "generated during a previous run ",
                            default=None)
        parser.add_argument("-v", "--verbose", dest="verbose",
                            action="store_true", default=False,
                            help="enable verbose logging output "
                            "[default=%(default)s]")
        parser.add_argument("-p", "--processors", dest="num_processors", 
                            type=int, default=DEFAULT_NUM_PROCESSORS,
                            help="Number of processor cores to allocate to "
                            "chimerascan [default=%(default)s]")
        parser.add_argument("--keep-tmp", dest="keep_tmp", 
                            action="store_true",
                            default=DEFAULT_KEEP_TMP,
                            help="DO NOT delete intermediate files after "
                            "run [default=%(default)s]")
        parser.add_argument("--rm-tmp", dest="keep_tmp", 
                            action="store_false", 
                            help="Delete intermediate files after run "
                            "[default=%s]" % str(not DEFAULT_KEEP_TMP))
        parser.add_argument("--quals", dest="quals",
                            choices=FASTQ_QUAL_FORMATS, 
                            default=DEFAULT_FASTQ_QUAL_FORMAT, metavar="FMT",
                            help="FASTQ quality score format "
                            "[default=%(default)s]")
        parser.add_argument('--library-type', dest="library_type", 
                            choices=LibraryTypes.choices(),
                            default=DEFAULT_LIBRARY_TYPE,
                            help="Library type [default=%(default)s]")
        parser.add_argument("--isize-mean", dest="isize_mean", type=int,
                            default=DEFAULT_ISIZE_MEAN, metavar="N",
                            help="Mean insert size to sample from when "
                            "insert size distribution cannot be determined "
                            "empirically [default=%(default)s]")
        parser.add_argument("--isize-stdev", dest="isize_stdev", type=float,
                            default=DEFAULT_ISIZE_STDEV, metavar="N",
                            help="Insert size standard deviation to sample "
                            "from when insert size distribution cannot be "
                            "determined empirically [default=%(default)s]")
        parser.add_argument("--trim5", type=int, dest="trim5", 
                            default=DEFAULT_TRIM5, metavar="N",
                            help="Trim N bases from 5' end of read")
        parser.add_argument("--trim3", type=int, dest="trim3", 
                            default=DEFAULT_TRIM3, metavar="N",
                            help="Trim N bases from 3' end of read")
        parser.add_argument("--min-fragment-length", type=int, 
                            dest="min_fragment_length", 
                            default=DEFAULT_MIN_FRAG_LENGTH,
                            help="Smallest expected fragment length "
                            "[default=%(default)s]")
        parser.add_argument("--max-fragment-length", type=int, 
                            dest="max_fragment_length", 
                            default=DEFAULT_MAX_FRAG_LENGTH,
                            help="Largest expected fragment length (reads "
                            "less than this fragment length are assumed to "
                            "be unspliced and contiguous) "
                            "[default=%(default)s]")
        parser.add_argument("--segment-length", type=int, 
                            dest="segment_length", 
                            default=DEFAULT_SEGMENT_LENGTH,
                            metavar="N",
                            help="Override size of soft-clipped read "
                            "segments during discordant alignment phase "
                            "(determined empirically by default)")
        parser.add_argument("--multihits", type=int, 
                            dest="max_multihits", 
                            default=DEFAULT_MAX_MULTIHITS,
                            metavar="N",
                            help="Override size of soft-clipped read "
                            "segments during discordant alignment phase "
                            "(determined empirically by default)")
        # filtering options
        group = parser.add_argument_group('Filtering options')
        group.add_argument("--homology-mismatches", type=int, 
                           dest="homology_mismatches",
                           default=DEFAULT_HOMOLOGY_MISMATCHES,
                           help="Number of mismatches to tolerate at "
                           "breakpoints when computing homology "
                           "[default=%(default)s]", metavar="N")
        group.add_argument("--anchor-min", type=int, dest="anchor_min", 
                           default=DEFAULT_ANCHOR_MIN, metavar="N",
                           help="Minimum breakpoint overlap (bp) required "
                           "to call spanning reads [default=%(default)s]")
        group.add_argument("--anchor-length", type=int, dest="anchor_length", 
                           default=DEFAULT_ANCHOR_LENGTH, metavar="N",
                           help="Size (bp) of anchor region where mismatch "
                           "checks are enforced [default=%(default)s]")
        group.add_argument("--anchor-mismatches", type=int, 
                           dest="anchor_mismatches", 
                           default=DEFAULT_ANCHOR_MISMATCHES,
                           metavar="N",
                           help="Number of mismatches allowed within anchor "
                           "region [default=%(default)s]")
        group.add_argument("--filter-unique-frags", type=float,
                           default=DEFAULT_FILTER_UNIQUE_FRAGS,
                           dest="filter_unique_frags", metavar="N",
                           help="Filter chimeras with less than N unique "
                           "aligned fragments [default=%(default)s]")
        group.add_argument("--filter-isize-prob", type=float,
                           default=DEFAULT_FILTER_ISIZE_PROB,
                           dest="filter_isize_prob", metavar="X",
                           help="Filter chimeras when probability of "
                           "observing the putative insert size is less "
                           "than X (0.0-1.0) [default=%(default)s]")
        group.add_argument("--filter-isoform-fraction", type=float, 
                           default=DEFAULT_FILTER_ISOFORM_FRACTION, metavar="X",
                           help="Filter chimeras with expression ratio "
                           "less than X (0.0-1.0) relative to the wild-type "
                           "transcript levels [default=%(default)s]")
        group.add_argument("--filter-false-pos", default="",
                           dest="filter_false_pos_file",
                           help="File containing known false positive "
                           "chimeric transcript pairs to filter out")
        return parser

    def from_args(self, args, parser=None):
        if parser is None:
            parser = self.get_argument_parser()
        args = parser.parse_args(args=args)
        # parse config file options/args
        if args.config_file is not None:
            self.from_xml(args.config_file)
        # reset logging to verbose
        if args.verbose:
            logging.getLogger().setLevel(logging.DEBUG)
        # check command line arguments
        if self.index_dir is None:
            if args.index_dir is None:
                parser.error("index not specified in config file or command line")
            self.index_dir = os.path.abspath(args.index_dir)
        if self.fastq_files[0] is None:
            self.fastq_files[0] = os.path.abspath(args.read1)
        if self.fastq_files[1] is None:
            self.fastq_files[1] = os.path.abspath(args.read2)
        if None in self.fastq_files:            
            parser.error("fastq files not specified in config file or command line")        
        if self.output_dir is None:
            if args.output_dir is None:
                parser.error("output dir not specified in config file or command line")   
            self.output_dir = os.path.abspath(args.output_dir)
        # optional arguments
        # set rest of options, overriding if attribute is undefined
        # or set to something other than the default 
        for attrname, attrtype, attrdefault in self.attrs:
            if ((getattr(self, attrname) is None) or
                (getattr(args, attrname) != attrdefault)):
                setattr(self, attrname, getattr(args, attrname))        

    def check_config(self):
        # check that input fastq files exist
        config_passed = True
        for mate,fastq_file in enumerate(self.fastq_files):
            if not os.path.isfile(fastq_file):
                logging.error("mate '%d' fastq file '%s' is not valid" % 
                              (mate, fastq_file))
                config_passed = False
        # check read lengths with trimming applied
        logging.debug("Checking read lengths")
        read_lengths = [detect_read_length(fq) for fq in self.fastq_files]
        total_trimming = self.trim5 + self.trim3
        for i,rlen in enumerate(read_lengths):
            trimmed_rlen = rlen - total_trimming
            logging.debug("File %s read length: %d after trimming: %d" % 
                          (self.fastq_files[i], rlen, trimmed_rlen))
            if trimmed_rlen < config.MIN_SEGMENT_LENGTH:
                logging.error("Trimmed read length is less than the minimum length of %d" % 
                              (trimmed_rlen, config.MIN_SEGMENT_LENGTH))
                config_passed = False
        # check that mate read lengths are equal
        if len(set(read_lengths)) > 1:
            logging.error("Unequal read lengths mate1=%d and mate2=%d" % 
                          (read_lengths[0], read_lengths[1]))
            config_passed = False
        # check that seed length < read length
        if self.segment_length is not None:
            if any((self.segment_length > rlen) for rlen in read_lengths):
                logging.error("seed length %d cannot be longer than read length" % 
                              (self.segment_length))
                config_passed = False
        # check that output dir is not a regular file
        if os.path.exists(self.output_dir) and (not os.path.isdir(self.output_dir)):
            logging.error("Output directory name '%s' exists and is not a valid directory" % 
                          (self.output_dir))
            config_passed = False
        if check_executable(config.BOWTIE2_BUILD_BIN):
            logging.debug("Checking for '%s' binary... found" % config.BOWTIE2_BUILD_BIN)
        else:
            logging.error("%s binary not found or not executable" % config.BOWTIE2_BUILD_BIN)
            config_passed = False
        # check that bowtie program exists
        if check_executable(os.path.join(config.BOWTIE2_BIN)):
            logging.debug("Checking for '%s' binary... found" % config.BOWTIE2_BIN)
        else:
            logging.error("%s binary not found or not executable" % config.BOWTIE2_BIN)
            config_passed = False
        # check that alignment index exists
        if os.path.isdir(self.index_dir):
            logging.debug("Checking for chimerascan index directory... found")
            # check that alignment index files exist
            for f in config.TRANSCRIPTOME_BOWTIE2_FILES:
                filename = os.path.join(self.index_dir, f)
                if not os.path.isfile(filename):
                    logging.error("chimerascan index file '%s' invalid" % (filename))
                    config_passed = False
                    break
            for f in config.GENOME_BOWTIE2_FILES:
                filename = os.path.join(self.index_dir, f)
                if not os.path.isfile(filename):
                    logging.error("chimerascan index file '%s' invalid" % (filename))
                    config_passed = False
                    break
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
        return config.JOB_ERROR
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
    #
    # Setup run parameters
    #
    transcript_file = os.path.join(runconfig.index_dir, config.TRANSCRIPT_FEATURE_FILE)
    genome_index = os.path.join(runconfig.index_dir, config.GENOME_INDEX)
    transcriptome_index = os.path.join(runconfig.index_dir, config.TRANSCRIPTOME_INDEX)
    maxhits_file = os.path.join(runconfig.index_dir, 
                                      config.MAX_MULTIMAPPING_FILE)
    maxhits = int(open(maxhits_file).next().strip())
    original_read_length = detect_read_length(runconfig.fastq_files[0])
    # minimum fragment length cannot be smaller than the trimmed read length
    trimmed_read_length = (original_read_length - runconfig.trim5 - runconfig.trim3)
    min_fragment_length = max(runconfig.min_fragment_length, trimmed_read_length)
    # 
    # Process and inspect the FASTQ files, performing several alterations 
    # to the reads:
    #
    # 1) rename them from long string to numbers to save space throughout
    #    the pipeline. also store mapping from read numbers to full names 
    #    in a separate file
    # 2) ensure the "/1" and "/2" suffixes exist to denote paired reads
    # 3) convert quality scores to sanger format
    # 
    converted_fastq_files = [os.path.join(tmp_dir, fq) 
                             for fq in config.CONVERTED_FASTQ_FILES]
    read_name_dbm_file = os.path.join(tmp_dir, config.READ_NAME_DBM_FILE)
    msg = "Processing FASTQ files"
    skip = all(up_to_date(cfq, fq) for cfq,fq in 
               zip(converted_fastq_files, runconfig.fastq_files))
    skip = skip and up_to_date(read_name_dbm_file, runconfig.fastq_files[0])
    if skip:
        logging.info("[SKIPPED] %s" % (msg))
    else:
        logging.info(msg)
        converted_fastq_prefix = \
            os.path.join(tmp_dir, config.CONVERTED_FASTQ_PREFIX)
        try:
            retcode = process_input_reads(runconfig.fastq_files, 
                                          converted_fastq_prefix,
                                          quals=runconfig.quals,
                                          trim5=runconfig.trim5,
                                          trim3=runconfig.trim3)
            if retcode != config.JOB_SUCCESS:
                logging.error("%s step failed" % (msg))
                return config.JOB_ERROR
        except Exception as e:
            logging.info("Cleaning up after error %s" % (str(e)))
            for fq in converted_fastq_files:
                if os.path.isfile(fq):
                    os.remove(fq)
    #
    # Transcriptome alignment step
    #
    # Align to transcriptome in paired-end mode, trying to resolve as many 
    # reads as possible.
    #
    transcriptome_bam_file = os.path.join(runconfig.output_dir, config.TRANSCRIPTOME_BAM_FILE)
    transcriptome_unaligned_path = os.path.join(tmp_dir, config.TRANSCRIPTOME_UNALIGNED_PATH)
    transcriptome_unaligned_fastq_files = tuple(os.path.join(tmp_dir, fq) for fq in config.TRANSCRIPTOME_UNALIGNED_FASTQ_FILES)
    msg = "Aligning paired-end reads to transcriptome"
    if (all(up_to_date(transcriptome_bam_file, fq) for fq in converted_fastq_files) and 
        all(up_to_date(a,b) for a,b in zip(transcriptome_unaligned_fastq_files, converted_fastq_files))):
        logging.info("[SKIPPED] %s" % (msg))
    else:
        logging.info(msg)
        log_file = os.path.join(log_dir, config.TRANSCRIPTOME_LOG_FILE)
        retcode = bowtie2_align_pe(index=transcriptome_index,
                                   fastq_files=converted_fastq_files,
                                   unaligned_path=transcriptome_unaligned_path,
                                   bam_file=transcriptome_bam_file,
                                   log_file=log_file,
                                   library_type=runconfig.library_type,
                                   min_fragment_length=min_fragment_length,
                                   max_fragment_length=runconfig.max_fragment_length,
                                   maxhits=maxhits,
                                   num_processors=runconfig.num_processors)
        # cleanup if job failed
        if retcode != config.JOB_SUCCESS:
            if os.path.exists(transcriptome_bam_file):
                os.remove(transcriptome_bam_file)
            for f in transcriptome_unaligned_fastq_files:
                if os.path.exists(f):
                    os.remove(f)
    #
    # Genome alignment step
    #
    # Align any unaligned transcriptome reads to genome in paired-end mode.
    # Resolve as many reads as possible.
    #
    genome_bam_file = os.path.join(runconfig.output_dir, config.GENOME_BAM_FILE)
    genome_unaligned_path = os.path.join(tmp_dir, config.GENOME_UNALIGNED_PATH)
    genome_unaligned_fastq_files = tuple(os.path.join(tmp_dir, fq) for fq in config.GENOME_UNALIGNED_FASTQ_FILES)
    msg = "Realigning unaligned paired-end reads to genome"
    if (all(up_to_date(genome_bam_file, fq) for fq in converted_fastq_files) and 
        all(up_to_date(a,b) for a,b in zip(genome_unaligned_fastq_files, converted_fastq_files))):
        logging.info("[SKIPPED] %s" % (msg))
    else:
        logging.info(msg)
        log_file = os.path.join(log_dir, config.GENOME_LOG_FILE)
        retcode = bowtie2_align_pe(index=genome_index,
                                   fastq_files=transcriptome_unaligned_fastq_files,
                                   unaligned_path=genome_unaligned_path,
                                   bam_file=genome_bam_file,
                                   log_file=log_file,
                                   library_type=runconfig.library_type,
                                   min_fragment_length=min_fragment_length,
                                   max_fragment_length=runconfig.max_fragment_length,
                                   maxhits=maxhits,
                                   num_processors=runconfig.num_processors)
        # cleanup if job failed
        if retcode != config.JOB_SUCCESS:
            if os.path.exists(genome_bam_file):
                os.remove(genome_bam_file)
            for f in genome_unaligned_fastq_files:
                if os.path.exists(f):
                    os.remove(f)
    # 
    # TODO: Convert transcriptome reads to genome, and index all genome 
    # reads together 
    # 
    #
    # Sort aligned reads by position
    #
    msg = "Sorting aligned reads"
    sorted_transcriptome_bam_file = os.path.join(runconfig.output_dir, 
                                                 config.SORTED_TRANSCRIPTOME_BAM_FILE)
    if (up_to_date(sorted_transcriptome_bam_file, transcriptome_bam_file)):
        logging.info("[SKIPPED] %s" % (msg))
    else:
        logging.info(msg)
        sorted_aligned_bam_prefix = os.path.splitext(sorted_transcriptome_bam_file)[0]
        pysam.sort("-m", str(int(1e9)), transcriptome_bam_file, sorted_aligned_bam_prefix)
    #
    # Index BAM file
    #
    msg = "Indexing BAM file"
    sorted_aligned_bam_index_file = sorted_transcriptome_bam_file + ".bai"
    if (up_to_date(sorted_aligned_bam_index_file, sorted_transcriptome_bam_file)):
        logging.info("[SKIPPED] %s" % (msg))
    else:
        logging.info(msg)
        pysam.index(sorted_transcriptome_bam_file)
    #
    # Get insert size distribution
    #
    isize_dist_file = os.path.join(runconfig.output_dir, 
                                   config.ISIZE_DIST_FILE)
    msg = "Profiling insert size distribution"
    if up_to_date(isize_dist_file, transcriptome_bam_file):
        logging.info("[SKIPPED] %s" % msg)
        isize_dist = InsertSizeDistribution.from_file(open(isize_dist_file, "r"))
    else:
        logging.info(msg)
        max_isize_samples = config.ISIZE_MAX_SAMPLES
        bamfh = pysam.Samfile(transcriptome_bam_file, "rb")
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
    #
    # Determine ideal segment length automatically
    #
    # log insert size statistics
    logging.info("Insert size samples=%d mean=%f std=%f median=%d mode=%d" % 
                 (isize_dist.n, isize_dist.mean(), isize_dist.std(), 
                  isize_dist.isize_at_percentile(50.0), isize_dist.mode()))    
    # choose a segment length to optimize mapping
    optimal_isize = isize_dist.isize_at_percentile(DEFAULT_FRAG_SIZE_SENSITIVITY)
    logging.info("Determining soft-clipped segment length")
    logging.debug("\tInsert size at %f percent of distribution is %d" % 
                 (DEFAULT_FRAG_SIZE_SENSITIVITY, optimal_isize))
    optimal_segment_length = int(round(optimal_isize / 2.0))
    logging.debug("\tOptimal segment length is %d" % (optimal_segment_length))
    segment_length = min(optimal_segment_length, trimmed_read_length)
    segment_length = max(config.MIN_SEGMENT_LENGTH, segment_length)
    logging.debug("\tAfter adjusting for min %d and read length %d, final segment length is %d" % 
                 (config.MIN_SEGMENT_LENGTH, trimmed_read_length, segment_length))
    if runconfig.segment_length is not None:
        logging.debug("\tOverriding auto segment length and using segment length of %d" % (runconfig.segment_length))
        segment_length = runconfig.segment_length
    #
    # Realignment step
    #
    # trim and realign all the initially unaligned reads in order to
    # increase sensitivity to detect reads spanning fusion junctions
    #
    realigned_bam_file = os.path.join(tmp_dir, config.REALIGNED_BAM_FILE)
    realigned_log_file = os.path.join(log_dir, config.REALIGNED_LOG_FILE)
    msg = "Trimming and realigning initially unmapped reads"
    if all(up_to_date(realigned_bam_file, fq) for fq in genome_unaligned_fastq_files):
        logging.info("[SKIPPED] %s" % (msg))
    else:
        logging.info(msg)           
        retcode = bowtie2_align_pe_sr(index=transcriptome_index,
                                      transcript_file=transcript_file,
                                      fastq_files=genome_unaligned_fastq_files,
                                      bam_file=realigned_bam_file,
                                      log_file=realigned_log_file,
                                      tmp_dir=tmp_dir,
                                      segment_length=segment_length,
                                      maxhits=maxhits,
                                      max_multihits=runconfig.max_multihits,
                                      num_processors=runconfig.num_processors)
        if retcode != config.JOB_SUCCESS:
            if os.path.exists(realigned_bam_file):
                os.remove(realigned_bam_file)
            return config.JOB_ERROR
    #
    # Find discordant reads
    #
    # iterate through realigned reads and divide them into groups of
    # concordant, discordant within a gene (isoforms), discordant
    # between different genes, and discordant in the genome
    #
    realigned_paired_bam_file = os.path.join(tmp_dir, config.REALIGNED_PAIRED_BAM_FILE)
    realigned_unmapped_bam_file = os.path.join(tmp_dir, config.REALIGNED_UNMAPPED_BAM_FILE)
    msg = "Classifying concordant and discordant read pairs"
    if (up_to_date(realigned_paired_bam_file, realigned_bam_file) and
        up_to_date(realigned_unmapped_bam_file, realigned_bam_file)):
        logging.info("[SKIPPED] %s" % (msg))
    else:
        logging.info(msg)
        retcode = find_discordant_fragments(realigned_bam_file, 
                                            paired_bam_file=realigned_paired_bam_file, 
                                            unmapped_bam_file=realigned_unmapped_bam_file,
                                            index_dir=runconfig.index_dir,
                                            max_isize=runconfig.max_fragment_length,
                                            library_type=runconfig.library_type)
        if retcode != config.JOB_SUCCESS:
            if os.path.exists(realigned_paired_bam_file):
                os.remove(realigned_paired_bam_file)
            if os.path.exists(realigned_unmapped_bam_file):
                os.remove(realigned_unmapped_bam_file)
            return config.JOB_ERROR
    #
    # Write a BEDPE file from discordant reads
    #
    discordant_bedpe_file = os.path.join(tmp_dir, config.DISCORDANT_BEDPE_FILE)
    msg = "Converting discordant reads to BEDPE format"
    if (up_to_date(discordant_bedpe_file, realigned_paired_bam_file)):
        logging.info("[SKIPPED] %s" % (msg))
    else:
        logging.info(msg)
        discordant_reads_to_bedpe(runconfig.index_dir,
                                  realigned_paired_bam_file,
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
                                    max_read_length=trimmed_read_length,
                                    homology_mismatches=runconfig.homology_mismatches)
        if retcode != config.JOB_SUCCESS:
            logging.error("Failed with error code %d" % (retcode))
            return config.JOB_ERROR
    #
    # Filter chimeras with few supporting fragments
    #
    filtered_encompassing_chimera_file = os.path.join(tmp_dir, config.FILTERED_ENCOMPASSING_CHIMERA_FILE)
    msg = "Filtering encompassing chimeras with few supporting reads"
    if (up_to_date(filtered_encompassing_chimera_file, encompassing_chimera_file)):
        logging.info("[SKIPPED] %s" % (msg))
    else:        
        logging.info(msg)
        retcode = filter_encompassing_chimeras(encompassing_chimera_file,
                                               filtered_encompassing_chimera_file,
                                               runconfig.filter_unique_frags)
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
    if (up_to_date(breakpoint_chimera_file, filtered_encompassing_chimera_file) and
        up_to_date(breakpoint_map_file, filtered_encompassing_chimera_file) and
        up_to_date(breakpoint_fasta_file, filtered_encompassing_chimera_file)):
        logging.info("[SKIPPED] %s" % (msg))
    else:
        logging.info(msg)
        chimeras_to_breakpoints(input_file=filtered_encompassing_chimera_file,
                                breakpoint_sorted_chimera_file=breakpoint_chimera_file, 
                                breakpoint_map_file=breakpoint_map_file, 
                                breakpoint_fasta_file=breakpoint_fasta_file,
                                tmp_dir=tmp_dir)
    #
    # Build a bowtie index to align and detect spanning reads
    #
    breakpoint_index = os.path.join(tmp_dir, config.BREAKPOINT_INDEX)
    breakpoint_index_files = (os.path.join(tmp_dir, f) for f in config.BREAKPOINT_BOWTIE2_FILES)
    msg = "Building bowtie index of breakpoint spanning sequences"
    skip = False
    if all(up_to_date(f, breakpoint_fasta_file) for f in breakpoint_index_files):
        logging.info("[SKIPPED] %s" % (msg))
    else:
        logging.info(msg)
        args = [config.BOWTIE2_BUILD_BIN, 
                breakpoint_fasta_file, 
                breakpoint_index]
        f = open(os.path.join(log_dir, "breakpoint_bowtie_index.log"), "w")
        retcode = subprocess.call(args, stdout=f, stderr=f)
        f.close()
        if retcode != config.JOB_SUCCESS:
            logging.error("'bowtie-build' failed to create breakpoint index")
            for f in breakpoint_index_files:
                if os.path.exists(f):
                    os.remove(f)
            return config.JOB_ERROR
    #
    # Extract reads that were encompassing when trimmed but may span
    # breakpoints when extended.
    #
    encomp_spanning_fastq_file = os.path.join(tmp_dir, config.ENCOMP_SPANNING_FASTQ_FILE)
    msg = "Extracting encompassing reads that may extend past breakpoints"
    if (up_to_date(encomp_spanning_fastq_file, filtered_encompassing_chimera_file)):
        logging.info("[SKIPPED] %s" % (msg))
    else:
        logging.info(msg)    
        retcode = nominate_encomp_spanning_reads(filtered_encompassing_chimera_file, 
                                                 encomp_spanning_fastq_file)
        if retcode != config.JOB_SUCCESS:
            logging.error("Failed to extract breakpoint-spanning reads")
            if os.path.exists(encomp_spanning_fastq_file):
                os.remove(encomp_spanning_fastq_file)
            return config.JOB_ERROR
    #
    # Process unmapped reads to predict spanning reads where one read
    # maps to a chimera and the other is unmapped
    #
    # First, annotate the single-mapped reads in the BAM file
    #
    single_mapped_bam_file = os.path.join(tmp_dir, config.SINGLE_MAPPED_BAM_FILE)
    unaligned_spanning_fastq_file = os.path.join(tmp_dir, config.UNALIGNED_SPANNING_FASTQ_FILE)
    msg = "Separating unmapped and single-mapping reads that may span breakpoints"
    if (up_to_date(single_mapped_bam_file, realigned_unmapped_bam_file) and
        up_to_date(single_mapped_bam_file, filtered_encompassing_chimera_file) and
        up_to_date(unaligned_spanning_fastq_file, filtered_encompassing_chimera_file, nzsize=False)):        
        logging.info("[SKIPPED] %s" % (msg))
    else:
        logging.info(msg)
        retcode = extract_single_mapped_reads(chimera_file=filtered_encompassing_chimera_file, 
                                              unmapped_bam_file=realigned_unmapped_bam_file,
                                              single_mapped_bam_file=single_mapped_bam_file,
                                              unmapped_fastq_file=unaligned_spanning_fastq_file,
                                              library_type=runconfig.library_type,
                                              tmp_dir=tmp_dir)
        if retcode != config.JOB_SUCCESS:
            logging.error("\tFailed")
            if os.path.exists(unaligned_spanning_fastq_file):
                os.remove(unaligned_spanning_fastq_file)
            if os.path.exists(single_mapped_bam_file):
                os.remove(single_mapped_bam_file)
            return config.JOB_ERROR
    #
    # Parse the single-mapped reads and choose putative spanning reads
    #
    msg = "Extracting single-mapped reads that may span breakpoints"
    single_mapped_spanning_fastq_file = os.path.join(tmp_dir, config.SINGLEMAP_SPANNING_FASTQ_FILE)
    if (up_to_date(single_mapped_spanning_fastq_file, filtered_encompassing_chimera_file, nzsize=False)):
        logging.info("[SKIPPED] %s" % (msg))
    else:
        logging.info(msg)
        retcode = nominate_single_mapped_spanning_reads(chimera_file=filtered_encompassing_chimera_file, 
                                                        single_mapped_bam_file=single_mapped_bam_file,
                                                        single_mapped_fastq_file=single_mapped_spanning_fastq_file,
                                                        tmp_dir=tmp_dir)
        if retcode != config.JOB_SUCCESS:
            logging.error("Failed to extract unmapped reads")
            if os.path.exists(single_mapped_spanning_fastq_file):
                os.remove(single_mapped_spanning_fastq_file)
            return config.JOB_ERROR
    #
    # Re-align encompassing reads that overlap breakpoint junctions
    # 
    encomp_spanning_bam_file = os.path.join(tmp_dir, config.ENCOMP_SPANNING_BAM_FILE)
    encomp_spanning_log_file = os.path.join(log_dir, "bowtie_encomp_spanning.log")
    msg = "Realigning encompassing reads to breakpoints"
    skip = all(up_to_date(encomp_spanning_bam_file, f) for f in breakpoint_index_files)
    skip = skip and up_to_date(encomp_spanning_bam_file, encomp_spanning_fastq_file)
    if skip:
        logging.info("[SKIPPED] %s" % (msg))
    else:
        logging.info(msg)
        retcode = bowtie2_align_sr(index=breakpoint_index,
                                   fastq_file=encomp_spanning_fastq_file,
                                   bam_file=encomp_spanning_bam_file,
                                   log_file=encomp_spanning_log_file,
                                   library_type=runconfig.library_type,
                                   maxhits=maxhits,
                                   num_processors=runconfig.num_processors)
        if retcode != config.JOB_SUCCESS:
            logging.error("Bowtie failed with error code %d" % (retcode))    
            if os.path.exists(encomp_spanning_bam_file):
                os.remove(encomp_spanning_bam_file)
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
    skip = all(up_to_date(singlemap_spanning_bam_file, f) for f in breakpoint_index_files)
    skip = skip and up_to_date(singlemap_spanning_bam_file, single_mapped_spanning_fastq_file)
    if skip:
        logging.info("[SKIPPED] %s" % (msg))
    else:            
        logging.info(msg)
        retcode = bowtie2_align_sr(index=breakpoint_index,
                                   fastq_file=single_mapped_spanning_fastq_file,
                                   bam_file=singlemap_spanning_bam_file,
                                   log_file=singlemap_spanning_log_file,
                                   library_type=runconfig.library_type,
                                   maxhits=maxhits,
                                   num_processors=runconfig.num_processors)
        if retcode != config.JOB_SUCCESS:
            logging.error("Bowtie failed with error code %d" % (retcode))
            if os.path.exists(singlemap_spanning_bam_file):
                os.remove(singlemap_spanning_bam_file)    
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
    # Merge spanning read alignment information
    #
    spanning_chimera_file = os.path.join(tmp_dir, config.SPANNING_CHIMERA_FILE)
    msg = "Merging spanning read information"
    if (up_to_date(spanning_chimera_file, breakpoint_chimera_file) and
        up_to_date(spanning_chimera_file, encomp_spanning_bam_file) and
        up_to_date(spanning_chimera_file, sorted_singlemap_spanning_bam_file)):
        logging.info("[SKIPPED] %s" % (msg))
    else:
        logging.info(msg)        
        merge_spanning_alignments(breakpoint_chimera_file=breakpoint_chimera_file,
                                  encomp_bam_file=sorted_encomp_spanning_bam_file,
                                  singlemap_bam_file=sorted_singlemap_spanning_bam_file,
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
    resolved_spanning_chimera_file = os.path.join(tmp_dir, config.RESOLVED_SPANNING_CHIMERA_FILE)
    msg = "Resolving ambiguous read mappings"
    if (up_to_date(resolved_spanning_chimera_file, spanning_chimera_file)):
        logging.info("[SKIPPED] %s" % (msg))
    else:
        logging.info(msg)        
        resolve_discordant_reads(input_file=spanning_chimera_file,
                                 output_file=resolved_spanning_chimera_file,
                                 isize_dist=isize_dist,
                                 min_isize_prob=runconfig.filter_isize_prob,
                                 tmp_dir=tmp_dir)
    #
    # Filter chimeras
    # 
    filtered_chimera_file = os.path.join(tmp_dir, config.FILTERED_CHIMERA_FILE)
    msg = "Filtering chimeras"
    if up_to_date(filtered_chimera_file, resolved_spanning_chimera_file):
        logging.info("[SKIPPED] %s" % (msg))
    else:
        logging.info(msg)
        # get insert size at prob    
        filter_chimeras(input_file=resolved_spanning_chimera_file, 
                        output_file=filtered_chimera_file,
                        index_dir=runconfig.index_dir,
                        bam_file=sorted_transcriptome_bam_file,
                        unique_frags=runconfig.filter_unique_frags,
                        isoform_fraction=runconfig.filter_isoform_fraction,
                        false_pos_file=runconfig.filter_false_pos_file)
    #
    # Filter homologous genes
    # 
    homolog_filtered_chimera_file = os.path.join(tmp_dir, config.HOMOLOG_FILTERED_CHIMERA_FILE)
    msg = "Filtering homologous chimeras"
    if up_to_date(homolog_filtered_chimera_file, filtered_chimera_file):
        logging.info("[SKIPPED] %s" % (msg))
    else:
        logging.info(msg)
        min_isize = isize_dist.isize_at_percentile(1.0)
        max_isize = isize_dist.isize_at_percentile(99.0)
        filter_homologous_genes(input_file=filtered_chimera_file, 
                                output_file=homolog_filtered_chimera_file,
                                index_dir=runconfig.index_dir,
                                homolog_segment_length=segment_length-1,
                                min_isize=min_isize,
                                max_isize=max_isize,
                                maxhits=maxhits,
                                num_processors=runconfig.num_processors,
                                tmp_dir=tmp_dir)
    #
    # Choose best isoform for chimeras that share the same breakpoint 
    # 
    best_isoform_chimera_file = os.path.join(tmp_dir, config.BEST_FILTERED_CHIMERA_FILE)
    msg = "Choosing best isoform for each chimera"
    if up_to_date(best_isoform_chimera_file, homolog_filtered_chimera_file):
        logging.info("[SKIPPED] %s" % (msg))
    else:
        logging.info(msg)
        retcode = filter_highest_coverage_isoforms(index_dir=runconfig.index_dir, 
                                                   input_file=homolog_filtered_chimera_file, 
                                                   output_file=best_isoform_chimera_file)
    #
    # Write user-friendly output file
    # 
    chimera_output_file = os.path.join(runconfig.output_dir, config.CHIMERA_OUTPUT_FILE)
    msg = "Writing chimeras to file %s" % (chimera_output_file)
    if up_to_date(chimera_output_file, best_isoform_chimera_file):
        logging.info("[SKIPPED] %s" % (msg))
    else:
        logging.info(msg)
        write_output(best_isoform_chimera_file,
                     bam_file=sorted_transcriptome_bam_file, 
                     output_file=chimera_output_file,
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
    return config.JOB_SUCCESS


def main():
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    # parse run parameters in config file and command line
    runconfig = RunConfig()
    runconfig.from_args(sys.argv[1:])
    # run chimerascan
    return run_chimerascan(runconfig)

if __name__ == '__main__':
    sys.exit(main())
