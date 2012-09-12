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

import sys
import logging

# check for python version 2.7.0 or greater
if sys.version_info < (2,7,1):
    sys.stderr.write("You need python 2.7.1 or later to run chimerascan\n")
    sys.exit(1)

# standard modules
import os
import subprocess
import shutil
import argparse
import xml.etree.ElementTree as etree

# third-party modules
import pysam

# local imports
import chimerascan.lib.config as config
from chimerascan.lib.base import LibraryTypes, check_executable, \
    parse_bool, indent_xml, up_to_date
from chimerascan.lib.seq import FASTQ_QUAL_FORMATS, SANGER_FORMAT, detect_read_length
from chimerascan.lib.fragment_size_distribution import InsertSizeDistribution
from chimerascan.lib.feature import TranscriptFeature

from chimerascan.pipeline.process_input_reads import process_input_reads
from chimerascan.pipeline.align_bowtie2 import bowtie2_align_transcriptome_pe, bowtie2_align_pe, bowtie2_align_pe_sr
from chimerascan.pipeline.find_discordant_reads import find_discordant_fragments
from chimerascan.pipeline.transcriptome_to_genome import transcriptome_to_genome
from chimerascan.pipeline.sam_to_bam import sam_to_bam
from chimerascan.pipeline.cluster_discordant_reads import cluster_discordant_reads
from chimerascan.pipeline.pair_clusters import pair_discordant_clusters
from chimerascan.pipeline.breakpoint_realignment import realign_across_breakpoints
from chimerascan.pipeline.process_spanning_alignments import process_spanning_alignments
from chimerascan.pipeline.filter_chimeras import filter_chimeras
from chimerascan.pipeline.write_output import write_output

# global default parameters
DEFAULT_NUM_PROCESSORS = config.BASE_PROCESSORS
DEFAULT_KEEP_TMP = True

# default sequencing data parameters
DEFAULT_FASTQ_QUAL_FORMAT = SANGER_FORMAT
DEFAULT_LIBRARY_TYPE = LibraryTypes.FR_UNSTRANDED
DEFAULT_ISIZE_MEAN = 200
DEFAULT_ISIZE_STDEV = 40
DEFAULT_FRAG_SIZE_SENSITIVITY = 5.0
# defaults for alignment
DEFAULT_TRIM5 = 0
DEFAULT_TRIM3 = 0
DEFAULT_SEGMENT_LENGTH = None

class RunConfig(object):
    
    attrs = (("num_processors", int, DEFAULT_NUM_PROCESSORS),
             ("keep_tmp", parse_bool, DEFAULT_KEEP_TMP),
             ("quals", str, DEFAULT_FASTQ_QUAL_FORMAT),
             ("library_type", str, DEFAULT_LIBRARY_TYPE),
             ("isize_mean", int, DEFAULT_ISIZE_MEAN),
             ("isize_stdev", float, DEFAULT_ISIZE_STDEV),
             ("min_fragment_length", int, config.DEFAULT_MIN_FRAG_LENGTH),
             ("max_fragment_length", int, config.DEFAULT_MAX_FRAG_LENGTH),
             ("trim5", int, DEFAULT_TRIM5),
             ("trim3", int, DEFAULT_TRIM3),
             ("segment_length", int, DEFAULT_SEGMENT_LENGTH),
             ("max_multihits", int, config.DEFAULT_MAX_MULTIHITS),
             ("local_multihits", int, config.DEFAULT_LOCAL_MULTIHITS),
             ("local_anchor_length", int, config.DEFAULT_LOCAL_ANCHOR_LENGTH),
             ("filter_num_frags", float, config.DEFAULT_FILTER_FRAGS),
             ("filter_allele_fraction", float, config.DEFAULT_FILTER_ALLELE_FRACTION),
             ("mask_biotypes_file", str, ""),
             ("mask_rnames_file", str, ""))

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
                            default=config.DEFAULT_MIN_FRAG_LENGTH,
                            help="Smallest expected fragment length "
                            "[default=%(default)s]")
        parser.add_argument("--max-fragment-length", type=int, 
                            dest="max_fragment_length", 
                            default=config.DEFAULT_MAX_FRAG_LENGTH,
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
                            default=config.DEFAULT_MAX_MULTIHITS,
                            metavar="N",
                            help="Maximum alignments allowed for each "
                            "discordant read")
        parser.add_argument("--local-multihits", type=int, 
                            dest="local_multihits", 
                            default=config.DEFAULT_LOCAL_MULTIHITS,
                            metavar="N",
                            help="Maximum alignments allowed for each "
                            "discordant read")
        parser.add_argument("--local-anchor-length", type=int, 
                            dest="local_anchor_length", 
                            default=config.DEFAULT_LOCAL_ANCHOR_LENGTH,
                            metavar="N",
                            help="Number of bases that read must span "
                            "on each side of a chimera to be considered "
                            "a valid breakpoint read")
        # filtering options
        group = parser.add_argument_group('Filtering options')
        group.add_argument("--filter-num-frags", type=float,
                           default=config.DEFAULT_FILTER_FRAGS,
                           dest="filter_num_frags", metavar="N",
                           help="Filter chimeras with less than N "
                           "aligned fragments [default=%(default)s]")
        group.add_argument("--filter-allele-fraction", type=float, 
                           default=config.DEFAULT_FILTER_ALLELE_FRACTION, 
                           dest="filter_allele_fraction", metavar="X",
                           help="Filter chimeras with expression less than "
                           "the specified fraction of the total expression "
                           "level [default=%(default)s")            
        group.add_argument("--mask-biotypes-file", default="",
                           dest="mask_biotypes_file",
                           help="File containing list of gene biotypes "
                           "to ignore (ex. pseudogenes, rRNA)")
        group.add_argument("--mask-rnames-file", default="",
                           dest="mask_rnames_file",
                           help="File containing list of reference names "
                           "to ignore (ex. MT or chrM)")
        # filtering options
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
        # ensure local anchor length is larger than minimum
        if self.local_anchor_length < config.LOCAL_ANCHOR_LENGTH_MIN:
            logging.error("Local anchor length of %d < %d" % 
                          (self.local_anchor_length, config.LOCAL_ANCHOR_LENGTH_MIN))
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
    # mask biotypes and references
    mask_biotypes = set()
    if runconfig.mask_biotypes_file:
        logging.info("Reading biotypes mask file")
        mask_biotypes.update([line.strip() for line in open(runconfig.mask_biotypes_file)])
        logging.info("\tread biotypes: %s" % (','.join(sorted(mask_biotypes))))
    mask_rnames = set()
    if runconfig.mask_rnames_file:
        logging.info("Reading references mask file")
        mask_rnames.update([line.strip() for line in open(runconfig.mask_rnames_file)])
        logging.info("\tread references: %s" % (','.join(sorted(mask_rnames))))
    # read transcripts
    logging.info("Reading transcript features")
    transcript_file = os.path.join(runconfig.index_dir, config.TRANSCRIPT_FEATURE_FILE)
    transcripts = list(TranscriptFeature.parse(open(transcript_file)))
    logging.info("\tread %d transcripts" % (len(transcripts)))
    # setup alignment indexes
    genome_index = os.path.join(runconfig.index_dir, config.GENOME_INDEX)
    transcriptome_index = os.path.join(runconfig.index_dir, config.TRANSCRIPTOME_INDEX)
    max_transcriptome_hits_file = os.path.join(runconfig.index_dir, 
                                               config.MAX_MULTIMAPPING_FILE)
    max_transcriptome_hits = int(open(max_transcriptome_hits_file).next().strip())
    # detect read length
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
    read_name_file = os.path.join(tmp_dir, config.READ_NAME_TXT_FILE)
    msg = "Processing FASTQ files"
    skip = all(up_to_date(cfq, fq) for cfq,fq in 
               zip(converted_fastq_files, runconfig.fastq_files))
    skip = skip and up_to_date(read_name_file, runconfig.fastq_files[0])
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
    transcriptome_bam_file = os.path.join(tmp_dir, config.TRANSCRIPTOME_BAM_FILE)
    transcriptome_unaligned_path = os.path.join(tmp_dir, config.TRANSCRIPTOME_UNALIGNED_PATH)
    transcriptome_unaligned_fastq_files = tuple(os.path.join(tmp_dir, fq) for fq in config.TRANSCRIPTOME_UNALIGNED_FASTQ_FILES)
    msg = "Aligning paired-end reads to transcriptome"
    if (all(up_to_date(transcriptome_bam_file, fq) for fq in converted_fastq_files) and 
        all(up_to_date(a,b) for a,b in zip(transcriptome_unaligned_fastq_files, converted_fastq_files))):
        logging.info("[SKIPPED] %s" % (msg))
    else:
        logging.info(msg)
        log_file = os.path.join(log_dir, config.TRANSCRIPTOME_LOG_FILE)
        retcode = bowtie2_align_transcriptome_pe(transcriptome_index=transcriptome_index,
                                                 genome_index=genome_index,
                                                 transcript_file=transcript_file,     
                                                 fastq_files=converted_fastq_files,
                                                 unaligned_path=transcriptome_unaligned_path,
                                                 bam_file=transcriptome_bam_file,
                                                 log_file=log_file,
                                                 library_type=runconfig.library_type,
                                                 min_fragment_length=min_fragment_length,
                                                 max_fragment_length=runconfig.max_fragment_length,
                                                 max_transcriptome_hits=max_transcriptome_hits,
                                                 num_processors=runconfig.num_processors)
        # cleanup if job failed
        if retcode != config.JOB_SUCCESS:
            logging.error("[FAILED] %s" % (msg))
            if os.path.exists(transcriptome_bam_file):
                os.remove(transcriptome_bam_file)
            for f in transcriptome_unaligned_fastq_files:
                if os.path.exists(f):
                    os.remove(f)
            return config.JOB_ERROR
    #
    # Sort transcriptome reads by position
    #
    msg = "Sorting transcriptome reads"
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
    sorted_transcriptome_bam_index_file = sorted_transcriptome_bam_file + ".bai"
    if (up_to_date(sorted_transcriptome_bam_index_file, sorted_transcriptome_bam_file)):
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
        bamfh = pysam.Samfile(sorted_transcriptome_bam_file, "rb")
        isize_dist = InsertSizeDistribution.from_genome_bam(bamfh, transcripts, 
                                                            min_isize=min_fragment_length, 
                                                            max_isize=runconfig.max_fragment_length, 
                                                            max_samples=config.ISIZE_MAX_SAMPLES)
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
                                                            samples=config.ISIZE_MAX_SAMPLES)
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
    optimal_segment_length = int(round(optimal_isize / 3.0))
    logging.debug("\tOptimal segment length is %d/3.0 = %d" % (optimal_isize, optimal_segment_length))
    segment_length = min(optimal_segment_length, trimmed_read_length)
    segment_length = max(config.MIN_SEGMENT_LENGTH, segment_length)
    logging.debug("\tAfter adjusting for min %d and read length %d, final segment length is %d" % 
                 (config.MIN_SEGMENT_LENGTH, trimmed_read_length, segment_length))
    if runconfig.segment_length is not None:
        logging.debug("\tOverriding auto segment length and using segment length of %d" % (runconfig.segment_length))
        segment_length = runconfig.segment_length
    #
    # Genome alignment step
    #
    # Align any unaligned transcriptome reads to genome in paired-end mode.
    # Resolve as many reads as possible.
    #
    genome_bam_file = os.path.join(tmp_dir, config.GENOME_BAM_FILE)
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
                                   max_hits=max_transcriptome_hits,
                                   num_processors=runconfig.num_processors)
        # cleanup if job failed
        if retcode != config.JOB_SUCCESS:
            logging.error("[FAILED] %s" % (msg))
            if os.path.exists(genome_bam_file):
                os.remove(genome_bam_file)
            for f in genome_unaligned_fastq_files:
                if os.path.exists(f):
                    os.remove(f)
            return config.JOB_ERROR
    #
    # Realignment step
    #
    # trim and realign all the initially unaligned reads in order to
    # increase sensitivity to detect reads spanning fusion junctions
    #
    realigned_bam_file = os.path.join(tmp_dir, config.REALIGNED_BAM_FILE)
    realigned_log_file = os.path.join(log_dir, config.REALIGNED_LOG_FILE)
    msg = "Trimming and realigning initially unmapped reads"
    if (all(up_to_date(realigned_bam_file, fq) for fq in genome_unaligned_fastq_files) and
        up_to_date(realigned_bam_file, isize_dist_file)):
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
                                      max_hits=max_transcriptome_hits,
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
    paired_bam_file = os.path.join(tmp_dir, config.PAIRED_BAM_FILE)
    discordant_bam_file = os.path.join(tmp_dir, config.DISCORDANT_BAM_FILE)
    unpaired_bam_file = os.path.join(tmp_dir, config.UNPAIRED_BAM_FILE)
    unmapped_bam_file = os.path.join(tmp_dir, config.UNMAPPED_BAM_FILE)
    multimap_bam_file = os.path.join(tmp_dir, config.MULTIMAP_BAM_FILE)
    unresolved_bam_file = os.path.join(tmp_dir, config.UNRESOLVED_BAM_FILE)
    output_files = (paired_bam_file, discordant_bam_file, unpaired_bam_file,
                    unmapped_bam_file, multimap_bam_file, unresolved_bam_file)
    msg = "Classifying concordant and discordant read pairs"
    if (all(up_to_date(f, realigned_bam_file) for f in output_files)):
        logging.info("[SKIPPED] %s" % (msg))
    else:
        logging.info(msg)
        retcode = find_discordant_fragments(transcripts=transcripts,
                                            input_bam_file=realigned_bam_file,
                                            paired_bam_file=paired_bam_file,
                                            discordant_bam_file=discordant_bam_file,
                                            unpaired_bam_file=unpaired_bam_file,
                                            unmapped_bam_file=unmapped_bam_file,
                                            multimap_bam_file=multimap_bam_file,
                                            unresolved_bam_file=unresolved_bam_file,
                                            max_isize=runconfig.max_fragment_length,
                                            max_multihits=runconfig.max_multihits,
                                            library_type=runconfig.library_type)
        if retcode != config.JOB_SUCCESS:
            logging.error("[FAILED] %s" % (msg))
            for f in output_files:
                if os.path.exists(f):
                    os.remove(f)
            return config.JOB_ERROR
    #
    # Convert discordant transcriptome reads to genome coordinates
    #
    discordant_genome_bam_file = os.path.join(tmp_dir, config.DISCORDANT_GENOME_BAM_FILE)
    msg = "Converting discordant transcriptome hits to genomic coordinates"
    if (up_to_date(discordant_genome_bam_file, discordant_bam_file)):
        logging.info("[SKIPPED] %s" % (msg))
    else:
        logging.info(msg)        
        discordant_genome_sam_file = os.path.join(tmp_dir, config.DISCORDANT_GENOME_SAM_FILE)
        retcode = transcriptome_to_genome(genome_index, transcripts, 
                                          input_file=discordant_bam_file, 
                                          output_file=discordant_genome_sam_file,
                                          library_type=runconfig.library_type,
                                          input_sam=False,
                                          output_sam=True)
        if retcode != config.JOB_SUCCESS:
            logging.error("[FAILED] %s" % (msg))
            if os.path.exists(discordant_genome_sam_file):
                os.remove(discordant_genome_sam_file)
            return config.JOB_ERROR
        retcode = sam_to_bam(discordant_genome_sam_file, discordant_genome_bam_file)
        if retcode != config.JOB_SUCCESS:
            logging.error("[FAILED] %s" % (msg))
            if os.path.exists(discordant_genome_bam_file):
                os.remove(discordant_genome_bam_file)
            return config.JOB_ERROR
        if os.path.exists(discordant_genome_sam_file):
            os.remove(discordant_genome_sam_file)
    #
    # Sort discordant reads by position
    #
    msg = "Sorting discordant BAM file"
    sorted_discordant_genome_bam_file = os.path.join(tmp_dir, config.SORTED_DISCORDANT_GENOME_BAM_FILE)
    if (up_to_date(sorted_discordant_genome_bam_file, discordant_genome_bam_file)):
        logging.info("[SKIPPED] %s" % (msg))
    else:
        logging.info(msg)
        bam_prefix = os.path.splitext(sorted_discordant_genome_bam_file)[0]
        pysam.sort("-m", str(int(1e9)), discordant_genome_bam_file, bam_prefix)
    #
    # Index BAM file
    #
    msg = "Indexing discordant BAM file"
    sorted_discordant_bam_index_file = sorted_discordant_genome_bam_file + ".bai"
    if (up_to_date(sorted_discordant_bam_index_file, sorted_discordant_genome_bam_file)):
        logging.info("[SKIPPED] %s" % (msg))
    else:
        logging.info(msg)
        pysam.index(sorted_discordant_genome_bam_file)
    #
    # Convert unpaired transcriptome reads to genome coordinates
    #
    unpaired_genome_bam_file = os.path.join(tmp_dir, config.UNPAIRED_GENOME_BAM_FILE)
    msg = "Converting unpaired transcriptome hits to genomic coordinates"
    if (up_to_date(unpaired_genome_bam_file, unpaired_bam_file)):
        logging.info("[SKIPPED] %s" % (msg))
    else:
        logging.info(msg)        
        unpaired_genome_sam_file = os.path.join(tmp_dir, config.UNPAIRED_GENOME_SAM_FILE)
        retcode = transcriptome_to_genome(genome_index, transcripts, 
                                          input_file=unpaired_bam_file, 
                                          output_file=unpaired_genome_sam_file,
                                          library_type=runconfig.library_type,
                                          input_sam=False,
                                          output_sam=True)
        if retcode != config.JOB_SUCCESS:
            logging.error("[FAILED] %s" % (msg))
            if os.path.exists(unpaired_genome_sam_file):
                os.remove(unpaired_genome_sam_file)
            return config.JOB_ERROR
        retcode = sam_to_bam(unpaired_genome_sam_file, unpaired_genome_bam_file)
        if retcode != config.JOB_SUCCESS:
            logging.error("[FAILED] %s" % (msg))
            if os.path.exists(unpaired_genome_bam_file):
                os.remove(unpaired_genome_bam_file)
            return config.JOB_ERROR
        if os.path.exists(unpaired_genome_sam_file):
            os.remove(unpaired_genome_sam_file)        
    #
    # Sort unpaired reads by position
    #
    msg = "Sorting unpaired BAM file"
    sorted_unpaired_genome_bam_file = os.path.join(tmp_dir, config.SORTED_UNPAIRED_GENOME_BAM_FILE)
    if (up_to_date(sorted_unpaired_genome_bam_file, unpaired_genome_bam_file)):
        logging.info("[SKIPPED] %s" % (msg))
    else:
        logging.info(msg)
        bam_prefix = os.path.splitext(sorted_unpaired_genome_bam_file)[0]
        pysam.sort("-m", str(int(1e9)), unpaired_genome_bam_file, bam_prefix)
    #
    # Index BAM file
    #
    msg = "Indexing unpaired BAM file"
    sorted_unpaired_bam_index_file = sorted_unpaired_genome_bam_file + ".bai"
    if (up_to_date(sorted_unpaired_bam_index_file, sorted_unpaired_genome_bam_file)):
        logging.info("[SKIPPED] %s" % (msg))
    else:
        logging.info(msg)
        pysam.index(sorted_unpaired_genome_bam_file)
    #
    # Cluster discordant reads into chimera candidates
    #
    cluster_file = os.path.join(tmp_dir, config.DISCORDANT_CLUSTER_FILE)
    cluster_shelve_file = \
        os.path.join(tmp_dir, config.DISCORDANT_CLUSTER_SHELVE_FILE)
    sorted_discordant_genome_cluster_bam_file = \
        os.path.join(runconfig.output_dir, 
                     config.SORTED_DISCORDANT_GENOME_CLUSTER_BAM_FILE)
    input_files = (sorted_discordant_genome_bam_file, 
                   sorted_unpaired_genome_bam_file)
    output_files = (cluster_file, cluster_shelve_file,                      
                    sorted_discordant_genome_cluster_bam_file)
    msg = "Clustering discordant reads"
    skip = True
    for input_file in input_files:
        for output_file in output_files:
            skip = skip and up_to_date(output_file, input_file)
    if skip:
        logging.info("[SKIPPED] %s" % (msg))
    else:
        logging.debug(msg)
        retcode = cluster_discordant_reads(discordant_bam_file=sorted_discordant_genome_bam_file, 
                                           unpaired_bam_file=sorted_unpaired_genome_bam_file, 
                                           concordant_bam_file=sorted_transcriptome_bam_file, 
                                           output_bam_file=sorted_discordant_genome_cluster_bam_file, 
                                           cluster_file=cluster_file,
                                           cluster_shelve_file=cluster_shelve_file)
        if retcode != config.JOB_SUCCESS:
            logging.error("[FAILED] %s" % (msg))
            for f in output_files:
                if os.path.exists(f):
                    os.remove(f)
    #
    # Pair discordant clusters
    #
    cluster_pair_file = \
        os.path.join(tmp_dir, config.DISCORDANT_CLUSTER_PAIR_FILE)
    msg = "Pairing discordant clusters"
    output_files = (cluster_pair_file,)
    if up_to_date(cluster_pair_file, sorted_discordant_genome_cluster_bam_file):
        logging.info("[SKIPPED] %s" % (msg))
    else:
        logging.debug(msg)
        retcode = pair_discordant_clusters(discordant_bam_file=sorted_discordant_genome_cluster_bam_file, 
                                           cluster_pair_file=cluster_pair_file, 
                                           tmp_dir=tmp_dir)
        if retcode != config.JOB_SUCCESS:
            logging.error("[FAILED] %s" % (msg))
            for f in output_files:
                if os.path.exists(f):
                    os.remove(f)
    #
    # Perform realignment across putative fusion breakpoints
    #
    breakpoint_bam_file = os.path.join(tmp_dir, config.BREAKPOINT_BAM_FILE)
    msg = "Realigning to find breakpoint-spanning reads"
    input_files = (sorted_discordant_genome_bam_file, 
                   sorted_unpaired_genome_bam_file, 
                   cluster_shelve_file, 
                   cluster_pair_file)
    output_files = (breakpoint_bam_file,)
    skip = True
    for inp in input_files:
        for outp in output_files:
            if not up_to_date(outp, inp):
                skip = False
    if skip:
        logging.info("[SKIPPED] %s" % (msg))
    else:
        logging.debug(msg)
        retcode = realign_across_breakpoints(index_dir=runconfig.index_dir,
                                             discordant_bam_file=sorted_discordant_genome_bam_file,
                                             unpaired_bam_file=sorted_unpaired_genome_bam_file,
                                             cluster_shelve_file=cluster_shelve_file,
                                             cluster_pair_file=cluster_pair_file,
                                             breakpoint_bam_file=breakpoint_bam_file,
                                             log_dir=log_dir,
                                             tmp_dir=tmp_dir,
                                             num_processors=runconfig.num_processors,
                                             local_anchor_length=runconfig.local_anchor_length,
                                             local_multihits=runconfig.local_multihits)
        if retcode != config.JOB_SUCCESS:
            logging.error("[FAILED] %s" % (msg))
            for f in output_files:
                if os.path.exists(f):
                    os.remove(f)
    #
    # Nominate breakpoint spanning reads (split reads)
    #
    spanning_sam_file = os.path.join(tmp_dir, config.SPANNING_SAM_FILE)
    spanning_bam_file = os.path.join(tmp_dir, config.SPANNING_BAM_FILE)
    spanning_cluster_pair_file = os.path.join(tmp_dir, config.SPANNING_CLUSTER_PAIR_FILE)
    msg = "Processing breakpoint-spanning alignments"
    input_files = (breakpoint_bam_file,
                   cluster_shelve_file, 
                   cluster_pair_file)
    output_files = (spanning_bam_file,
                    spanning_cluster_pair_file)
    skip = True
    for inp in input_files:
        for outp in output_files:
            if not up_to_date(outp, inp):
                skip = False
    if skip:
        logging.info("[SKIPPED] %s" % (msg))
    else:
        logging.info(msg)
        retcode = process_spanning_alignments(cluster_shelve_file=cluster_shelve_file,
                                              cluster_pair_file=cluster_pair_file,
                                              bam_file=breakpoint_bam_file,                                              
                                              output_sam_file=spanning_sam_file,
                                              output_cluster_pair_file=spanning_cluster_pair_file,
                                              local_anchor_length=runconfig.local_anchor_length)
        if retcode != config.JOB_SUCCESS:
            logging.error("[FAILED] %s" % (msg))
            for f in output_files:
                if os.path.exists(f):
                    os.remove(f)
        retcode = sam_to_bam(spanning_sam_file, spanning_bam_file)
        if retcode != config.JOB_SUCCESS:
            logging.error("[FAILED] %s" % (msg))
            if os.path.exists(spanning_bam_file):
                os.remove(spanning_bam_file)
            return config.JOB_ERROR
        if os.path.exists(spanning_sam_file):
            os.remove(spanning_sam_file)
    #
    # Sort unpaired reads by position
    #
    msg = "Sorting spanning BAM file"
    sorted_spanning_bam_file = os.path.join(runconfig.output_dir, config.SORTED_SPANNING_BAM_FILE)
    if (up_to_date(sorted_spanning_bam_file, spanning_bam_file)):
        logging.info("[SKIPPED] %s" % (msg))
    else:
        logging.info(msg)
        bam_prefix = os.path.splitext(sorted_spanning_bam_file)[0]
        pysam.sort("-m", str(int(1e9)), spanning_bam_file, bam_prefix)
    #
    # Index BAM file
    #
    msg = "Indexing spanning BAM file"
    sorted_spanning_bam_index_file = sorted_spanning_bam_file + ".bai"
    if (up_to_date(sorted_spanning_bam_index_file, sorted_spanning_bam_file)):
        logging.info("[SKIPPED] %s" % (msg))
    else:
        logging.info(msg)
        pysam.index(sorted_spanning_bam_file)
    #
    # Write chimera file
    # 
    unfiltered_chimera_bedpe_file = os.path.join(runconfig.output_dir, 
                                                 config.UNFILTERED_CHIMERA_BEDPE_FILE)
    msg = "Writing unfiltered chimeras to file %s" % (unfiltered_chimera_bedpe_file)
    if (up_to_date(unfiltered_chimera_bedpe_file, spanning_cluster_pair_file) and
        up_to_date(unfiltered_chimera_bedpe_file, cluster_shelve_file)):                
        logging.info("[SKIPPED] %s" % (msg))
    else:
        logging.info(msg)
        retcode = write_output(transcripts, 
                               cluster_shelve_file=cluster_shelve_file, 
                               cluster_pair_file=spanning_cluster_pair_file, 
                               read_name_file=read_name_file, 
                               output_file=unfiltered_chimera_bedpe_file, 
                               annotation_source="ensembl")
        if retcode != config.JOB_SUCCESS:
            logging.error("[FAILED] %s" % (msg))
            if os.path.exists(unfiltered_chimera_bedpe_file):
                os.remove(unfiltered_chimera_bedpe_file)
    #
    # Filter chimeras
    #
    chimera_bedpe_file = os.path.join(runconfig.output_dir, config.CHIMERA_BEDPE_FILE)
    msg = "Filtering chimeras"
    if (up_to_date(chimera_bedpe_file, unfiltered_chimera_bedpe_file)):
        logging.info("[SKIPPED] %s" % (msg))
    else:        
        logging.info(msg)
        retcode = filter_chimeras(input_file=unfiltered_chimera_bedpe_file, 
                                  output_file=chimera_bedpe_file,
                                  filter_num_frags=runconfig.filter_num_frags,
                                  filter_allele_fraction=runconfig.filter_allele_fraction,
                                  mask_biotypes=mask_biotypes,
                                  mask_rnames=mask_rnames)
        if retcode != config.JOB_SUCCESS:
            logging.error("[FAILED] %s" % (msg))
            if os.path.exists(chimera_bedpe_file):
                os.remove(chimera_bedpe_file)
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
