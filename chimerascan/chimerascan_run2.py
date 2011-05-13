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
__author__ = "Matthew Iyer"
__copyright__ = "Copyright 2011, chimerascan project"
__credits__ = ["Matthew Iyer", "Christopher Maher"]
__license__ = "GPL"
__version__ = "0.3.4"
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
import chimerascan.pysam as pysam
import chimerascan.lib.config as config
from chimerascan.lib.config import JOB_SUCCESS, JOB_ERROR, indent
from chimerascan.lib.base import up_to_date, check_executable, \
    get_read_length, parse_bool
from chimerascan.lib.fragment_size_distribution import FragmentSizeDistribution

from chimerascan.pipeline.fastq_sort import sort_fastq_files
from chimerascan.pipeline.align_bowtie import align_pe_full
from chimerascan.pipeline.align_tophat import align_tophat


#from chimerascan.pipeline.align_full import align_pe_full, align_sr_full
#from chimerascan.pipeline.align_segments import align, determine_read_segments
#from chimerascan.pipeline.merge_read_pairs import merge_read_pairs
#from chimerascan.pipeline.find_discordant_reads import find_discordant_reads
#from chimerascan.pipeline.extend_sequences import extend_sequences
#from chimerascan.pipeline.sort_discordant_reads import sort_discordant_reads
#from chimerascan.pipeline.nominate_chimeras import nominate_chimeras
#from chimerascan.pipeline.filter_chimeras import filter_encompassing_chimeras
#from chimerascan.pipeline.bedpe_to_fasta import bedpe_to_junction_fasta
#from chimerascan.pipeline.merge_spanning_alignments import merge_spanning_alignments
#from chimerascan.pipeline.profile_insert_size import InsertSizeDistribution
#from chimerascan.pipeline.filter_spanning_chimeras import filter_spanning_chimeras
#from chimerascan.pipeline.rank_chimeras import rank_chimeras

# run defaults
DEFAULT_NUM_PROCESSORS = config.BASE_PROCESSORS
DEFAULT_KEEP_TMP = False

# alignment defaults
DEFAULT_TOPHAT_DIR = ""
DEFAULT_BOWTIE_DIR = ""
DEFAULT_SAMTOOLS_BIN = "samtools"
DEFAULT_INSERT_SIZE = 200
DEFAULT_INSERT_SIZE_STDEV = 20
DEFAULT_SEGMENT_LENGTH = 25
DEFAULT_MULTIHITS = 100
DEFAULT_MISMATCHES = 2
DEFAULT_TRIM5 = 0
DEFAULT_TRIM3 = 0
DEFAULT_MIN_FRAG_LENGTH = 0
DEFAULT_MAX_FRAG_LENGTH = 1000 
DEFAULT_MAX_INDEL_SIZE = 100
DEFAULT_QUALS = "sanger"
DEFAULT_LIBRARY_TYPE = "fr-unstranded"

DEFAULT_FILTER_MAX_MULTIMAPS = 2
DEFAULT_FILTER_MULTIMAP_RATIO = 0.10
DEFAULT_FILTER_STRAND_PVALUE = 0.01
DEFAULT_FILTER_ISIZE_PERCENTILE = 99.9

DEFAULT_ANCHOR_MIN = 0
DEFAULT_ANCHOR_MAX = 5
DEFAULT_ANCHOR_MISMATCHES = 5

DEFAULT_EMPIRICAL_PROB = 1.0

translate_quals = {'solexa': 'solexa-quals',
                   'illumina': 'solexa1.3-quals',
                   'sanger': 'phred33-quals'}
rev_translate_quals = dict((v,k) for k,v in translate_quals.items())

quals_help_string = ",".join(translate_quals.keys())

class RunConfig(object):
    
    attrs = (("num_processors", int, DEFAULT_NUM_PROCESSORS),
             ("keep_tmp", parse_bool, DEFAULT_KEEP_TMP),
             ("tophat_dir", str, DEFAULT_TOPHAT_DIR),
             ("bowtie_dir", str, DEFAULT_BOWTIE_DIR),
             ("samtools_bin", str, DEFAULT_SAMTOOLS_BIN),
             ("insert_size", int, DEFAULT_INSERT_SIZE),
             ("insert_size_stdev", int, DEFAULT_INSERT_SIZE_STDEV),
             ("segment_length", int, DEFAULT_SEGMENT_LENGTH),
             ("multihits", int, DEFAULT_MULTIHITS),
             ("mismatches", int, DEFAULT_MISMATCHES),
             ("trim5", int, DEFAULT_TRIM5),
             ("trim3", int, DEFAULT_TRIM3),
             ("min_fragment_length", int, DEFAULT_MIN_FRAG_LENGTH),
             ("max_fragment_length", int, DEFAULT_MAX_FRAG_LENGTH),
             ("max_indel_size", int, DEFAULT_MAX_INDEL_SIZE),
             ("library_type", str, DEFAULT_LIBRARY_TYPE),
             ("filter_max_multimaps", int, DEFAULT_FILTER_MAX_MULTIMAPS),
             ("filter_multimap_ratio", float, DEFAULT_FILTER_MULTIMAP_RATIO),
             ("filter_isize_percentile", float, DEFAULT_FILTER_ISIZE_PERCENTILE),
             ("filter_strand_pval", float, DEFAULT_FILTER_STRAND_PVALUE),
             ("anchor_min", int, DEFAULT_ANCHOR_MIN),
             ("anchor_max", int, DEFAULT_ANCHOR_MAX),
             ("anchor_mismatches", int, DEFAULT_ANCHOR_MISMATCHES),
             ("empirical_prob", float, DEFAULT_EMPIRICAL_PROB))

    def __init__(self):
        self.output_dir = None
        self.fastq_files = None
        self.index_dir = None
        self.quals = None
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
        quals = root.findtext("quals")
        if quals not in translate_quals:
            logging.error("Quality score option %s unknown, using default %s" % 
                          (quals, DEFAULT_QUALS))
        self.quals = quals
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
        # fastq format
        elem = etree.SubElement(root, "quals")
        elem.text = self.quals
        # optional arguments
        for attrname, attrtype, attrdefault in self.attrs:
            val = getattr(self, attrname)
            elem = etree.SubElement(root, attrname)
            elem.text = str(val)
        # indent for pretty printing
        indent(root)        
        return etree.tostring(root)

    @staticmethod
    def get_option_parser():
        parser = OptionParser(usage="%prog [options] [--config <config_file> "
                              " | <index_dir> <read1.fq> <read2.fq> "
                              "<output_dir>]",
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
                          help="Do not delete intermediate files from run")
        # alignment options
        align_group = OptionGroup(parser, "Alignment options",
                                  "Adjust these options to change " 
                                  "alignment settings")
        align_group.add_option("--tophat-dir", dest="tophat_dir", 
                               default=DEFAULT_TOPHAT_DIR, 
                               help="Path to directory containing 'tophat' " 
                               "binaries")
        align_group.add_option("--bowtie-dir", dest="bowtie_dir", 
                               default=DEFAULT_BOWTIE_DIR, 
                               help="Path to directory containing 'bowtie' " 
                               "binaries")
        align_group.add_option("--samtools-bin", dest="samtools_bin", 
                               default=DEFAULT_SAMTOOLS_BIN, 
                               help="Path to the 'samtools' binary")
        align_group.add_option("--multihits", type="int", dest="multihits", 
                               default=DEFAULT_MULTIHITS, metavar="MMAX",
                               help="Ignore reads that map to more than MMAX "
                               "locations [default=%default]")
        align_group.add_option("--mismatches", type="int", dest="mismatches",
                               default=DEFAULT_MISMATCHES, metavar="N",
                               help="Aligned reads must have <= N mismatches "
                               "[default=%default]")
        align_group.add_option("--segment-length", type="int", dest="segment_length", 
                               default=DEFAULT_SEGMENT_LENGTH,
                               help="Size of read segments (recommend 25 " 
                               "or higher for high specificity) "
                               "[default=%default]")
        align_group.add_option("--insert-size", type="int", dest="insert_size", 
                               default=DEFAULT_INSERT_SIZE,
                               help="Estimate total read fragment size "
                               "(beginning of read1 to end of read2) "
                               "[default=%default]")
        align_group.add_option("--insert-size-stdev", type="int", 
                               dest="insert_size_stdev", 
                               default=DEFAULT_INSERT_SIZE_STDEV,
                               help="Standard deviation of insert size "
                               "[default=%default]")
        align_group.add_option("--trim5", type="int", dest="trim5", 
                               default=DEFAULT_TRIM5, metavar="N",
                               help="Trim N bases from 5' end of read")
        align_group.add_option("--trim3", type="int", dest="trim3", 
                               default=DEFAULT_TRIM3, metavar="N",
                               help="Trim N bases from 3' end of read")
        align_group.add_option("--quals", dest="quals", 
                               default=DEFAULT_QUALS, metavar="FMT",
                               help="Choose from %s [default=%s]" % 
                               (quals_help_string, DEFAULT_QUALS))
        align_group.add_option("--min-fragment-length", type="int", 
                               dest="min_fragment_length", 
                               default=DEFAULT_MIN_FRAG_LENGTH,
                               help="Smallest expected fragment length "
                               "[default=%default]")
        align_group.add_option("--max-fragment-length", type="int", 
                               dest="max_fragment_length", 
                               default=DEFAULT_MAX_FRAG_LENGTH,
                               help="Largest expected fragment length (reads less"
                               " than this fragment length are assumed to be "
                               " unspliced and contiguous) [default=%default]")
        align_group.add_option('--library-type', dest="library_type",
                               choices=("fr-unstranded", 
                                        "fr-firststrand", 
                                        "fr-secondstrand"), 
                               default=DEFAULT_LIBRARY_TYPE,
                               help="Library type (see Tophat manual for "
                               " more information[default=%default]")        
        parser.add_option_group(align_group)
        # filtering options
        filter_group = OptionGroup(parser, "Filtering options",
                                   "Adjust these options to change "
                                   "filtering behavior")
        filter_group.add_option("--max-indel-size", type="int", 
                                dest="max_indel_size", 
                                default=DEFAULT_MAX_INDEL_SIZE,
                                help="Tolerate indels less than N bp "
                                "[default=%default]", metavar="N")
        filter_group.add_option("--anchor-min", type="int", dest="anchor_min", 
                                default=DEFAULT_ANCHOR_MIN,
                                help="Minimum junction overlap required to call "
                                "spanning reads [default=%default]")
        filter_group.add_option("--anchor-max", type="int", dest="anchor_max", 
                                default=DEFAULT_ANCHOR_MAX,
                                help="Junction overlap below which to enforce "
                                "mismatch checks [default=%default]")
        filter_group.add_option("--anchor-mismatches", type="int", dest="anchor_mismatches", 
                                default=DEFAULT_ANCHOR_MISMATCHES,
                                help="Number of mismatches allowed within anchor "
                                "region [default=%default]")        
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
        filter_group.add_option("--filter-isize-prob", type="int",
                                default=DEFAULT_FILTER_ISIZE_PERCENTILE,
                                dest="filter_isize_percentile", metavar="N",
                                help="Filter chimeras when putative insert "
                                "size is larger than the (1-N) percentile "
                                "of the distribution [default=%default]")
        filter_group.add_option("--filter-strand-pvalue", type="float",
                                default=DEFAULT_FILTER_STRAND_PVALUE,
                                dest="filter_strand_pval", metavar="p",
                                help="p-value to reject chimera based on "
                                " binomial test that balance of +/- strand "
                                " encompassing reads should be 50/50 "
                                "[default=%default]")
        filter_group.add_option("--empirical-prob", type="float", metavar="p", 
                                dest="empirical_prob", default=DEFAULT_EMPIRICAL_PROB, 
                                help="empirical probability threshold for "
                                "outputting chimeras [default=%default]")        
        parser.add_option_group(filter_group)
        return parser        

    def from_args(self, args, parser=None):
        if parser is None:
            parser = self.get_option_parser()
        options, args = parser.parse_args(args=args)
        # parse config file options/args
        if options.config_file is not None:
            self.from_xml(options.config_file)
        elif len(args) < 4:
            parser.error("expecting 'index_dir', 'read1', 'read2', and "
                         "'output_dir' on command line, but got less than "
                         "4 positional arguments")             
        # reset logging to verbose
        if options.verbose:
            logging.getLogger().setLevel(logging.DEBUG)
        # check command line arguments
        if self.index_dir is None:
            self.index_dir = args[0] 
        if self.fastq_files is None:
            self.fastq_files = args[1:3]
        if self.output_dir is None:
            self.output_dir = args[3]
        if self.quals is None:
            self.quals = options.quals
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
        # check that segment length < read length
        if any(self.segment_length > rlen for rlen in read_lengths):
            logging.error("segment length %d cannot be longer than read length" % 
                         (self.segment_length))
            config_passed = False
        # check that output dir is not a regular file
        if os.path.exists(self.output_dir) and (not os.path.isdir(self.output_dir)):
            logging.error("Output directory name '%s' exists and is not a valid directory" % 
                          (self.output_dir))
            config_passed = False
        if check_executable(os.path.join(self.bowtie_dir, "bowtie-build")):
            logging.debug("Checking for 'bowtie-build' binary... found")
        else:
            logging.error("'bowtie-build' binary not found or not executable")
            config_passed = False
        # check that bowtie program exists
        if check_executable(os.path.join(self.bowtie_dir, "bowtie")):
            logging.debug("Checking for 'bowtie' binary... found")
        else:
            logging.error("'bowtie' binary not found or not executable")
            config_passed = False
        # check the the tophat program exists
        if check_executable(os.path.join(self.tophat_dir, "tophat")):
            logging.debug("Checking for 'tophat' binary... found")
        else:
            logging.error("'tophat' binary not found or not executable")
            config_passed = False
        # check the the samtools program exists
        if check_executable(self.samtools_bin):
            logging.debug("Checking for 'samtools' binary... found")
        else:
            logging.error("'samtools' binary not found or not executable")
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
    # write the run config to a file
    xmlstring = runconfig.to_xml()
    runconfig_xml_file = os.path.join(runconfig.output_dir, config.RUNCONFIG_XML_FILE)
    fh = open(runconfig_xml_file, "w")
    print >>fh, xmlstring
    fh.close()
    # gather and parse run parameters
    gene_feature_file = os.path.join(runconfig.index_dir, config.GENE_FEATURE_FILE)
    bowtie_index = os.path.join(runconfig.index_dir, config.ALIGN_INDEX)
    bowtie_frag_size_index = os.path.join(runconfig.index_dir, config.FRAG_SIZE_INDEX)    
    tophat_juncs_file = os.path.join(runconfig.index_dir, config.TOPHAT_JUNCS_FILE)
    original_read_length = get_read_length(runconfig.fastq_files[0])
    # minimum fragment length cannot be smaller than the trimmed read length
    trimmed_read_length = original_read_length - runconfig.trim5 - runconfig.trim3
    min_fragment_length = max(runconfig.min_fragment_length, 
                              trimmed_read_length)
    #
    # Read processing step
    #
    # TODO: trim reads based on quality, etc


    #
    # Convert fastq files to BAM and sort
    #
    sorted_fastq_files = [os.path.join(tmp_dir, fq) 
                          for fq in config.SORTED_FASTQ_FILES]
    if all(up_to_date(sfq, fq) for sfq,fq in 
           zip(sorted_fastq_files, runconfig.fastq_files)):
        logging.info("[SKIPPED] Sorting fastq files by read name and "
                     "converting quality scores to 'sanger' format")
    else:
        logging.info("Sorting fastq files by read name and "
                     "converting quality scores to 'sanger' format")
        # convert fastq to BAM
        retcode = sort_fastq_files(runconfig.fastq_files, 
                                   sorted_fastq_files, 
                                   runconfig.quals, 
                                   tmp_dir)
        if retcode != 0:
            logging.error("Sort fastq files failed with error code %d" % (retcode))    
            return retcode

    #    
    # Align reads to special index to find fragment size distribution
    #
    frag_size_bam_file = os.path.join(tmp_dir, config.FRAG_SIZE_BAM_FILE)    
    frag_size_log_file = os.path.join(log_dir, "bowtie_frag_size_alignment.log")    
    msg = "Aligning reads to determine fragment size distribution"
    if all(up_to_date(frag_size_bam_file, fq) for fq in sorted_fastq_files):
        logging.info("[SKIPPED] %s" % (msg))
    else:
        logging.info(msg)
        # align to fragment size index
        bowtie_bin = os.path.join(runconfig.bowtie_dir, "bowtie")        
        retcode = align_pe_full(sorted_fastq_files, 
                                bowtie_frag_size_index,
                                frag_size_bam_file,
                                unaligned_fastq_param=None,
                                maxmultimap_fastq_param=None,
                                min_fragment_length=min_fragment_length,
                                max_fragment_length=runconfig.max_fragment_length,
                                trim5=runconfig.trim5,
                                trim3=runconfig.trim3,
                                library_type=runconfig.library_type,
                                num_processors=runconfig.num_processors,
                                quals="sanger",
                                multihits=1,
                                mismatches=2,
                                bowtie_bin=bowtie_bin,
                                bowtie_mode="-v",
                                log_file=frag_size_log_file)
        if retcode != 0:
            logging.error("Bowtie failed with error code %d" % (retcode))    
            sys.exit(retcode)
    #    
    # Estimate fragment size distribution
    #
    # this step uses the unambiguously mapping reads to estimate the
    # distribution of insert sizes in the library.  this will play a 
    # role in selecting an optimal trimmed read length to use when
    # remapping initially unmapped reads, and will also help to select
    # the most probable discordant read pairs when there are multiple 
    # mappings
    #
    frag_size_dist_file = os.path.join(runconfig.output_dir, config.FRAG_SIZE_DIST_FILE)
    if all(up_to_date(frag_size_dist_file, fq) for fq in sorted_fastq_files):
        logging.info("[SKIPPED] %s" % (msg))
        frag_size_dist = FragmentSizeDistribution.from_file(open(frag_size_dist_file, "r"))
    else:
        bamfh = pysam.Samfile(frag_size_bam_file, "rb")
        frag_size_dist = FragmentSizeDistribution.from_bam(bamfh, 
                            min_isize=min_fragment_length, 
                            max_isize=runconfig.max_fragment_length, 
                            max_samples=config.FRAG_SIZE_MAX_SAMPLES)
        frag_size_dist.to_file(open(frag_size_dist_file, "w"))
        bamfh.close()
    logging.info("Fragment size samples=%d mean=%f std=%f median=%d mode=%d" % 
                 (frag_size_dist.n, frag_size_dist.mean(), frag_size_dist.std(), 
                  frag_size_dist.percentile(50.0), frag_size_dist.mode()))
    #
    # Initial alignment step
    #
    # align in paired-end mode using Tophat, trying to resolve as many reads 
    # as possible this effectively rules out the vast majority of reads as 
    # candidate fusions and accounts for potential novel splice junctions as
    # well
    tophat_out_dir = os.path.join(runconfig.output_dir, config.TOPHAT_DIR)
    tophat_bam_file = os.path.join(tophat_out_dir, "accepted_hits.bam")
    tophat_log_file = os.path.join(log_dir, "tophat_alignment.log")
    msg = "Aligning reads in paired-end mode"
    if all(up_to_date(tophat_bam_file, fq) for fq in sorted_fastq_files):
        logging.info("[SKIPPED] %s" % (msg))
    else:    
        logging.info(msg)
        tophat_bin = os.path.join(runconfig.tophat_dir, "tophat")        
        retcode = align_tophat(sorted_fastq_files, 
                               tophat_bin, bowtie_index, tophat_out_dir,
                               read_length=trimmed_read_length,
                               insert_size=runconfig.insert_size,
                               insert_size_stdev=runconfig.insert_size_stdev,
                               num_processors=runconfig.num_processors,
                               library_type=runconfig.library_type,
                               quals="sanger",
                               multihits=runconfig.multihits,
                               mismatches=runconfig.mismatches,
                               segment_length=runconfig.segment_length,
                               raw_juncs_file=tophat_juncs_file,
                               log_file=tophat_log_file)
        if retcode != 0:
            logging.error("Bowtie failed with error code %d" % (retcode))    
            return retcode
    #    
    # Sort aligned reads by read ID 
    #
    # this step will sort the aligned reads by name, allowing iteration 
    # through the paired alignments to determine discordant reads
    #
    sorted_tophat_bam_file = os.path.join(tmp_dir, "rname_sorted_hits.bam")
    if up_to_date(sorted_tophat_bam_file, tophat_bam_file):
        logging.info("[SKIPPED] Sorting Tophat alignments by read name")
    else:
        logging.info("Sorting Tophat alignments by read name")
        pysam.sort("-n", tophat_bam_file, 
                   os.path.splitext(sorted_tophat_bam_file)[0])
    return
    



    
#    #
#    # Initial Bowtie alignment step
#    #
#    # align in paired-end mode, trying to resolve as many reads as possible
#    # this effectively rules out the vast majority of reads as candidate
#    # fusions
#    unaligned_fastq_param = os.path.join(tmp_dir, config.UNALIGNED_FASTQ_PARAM)
#    maxmultimap_fastq_param = os.path.join(tmp_dir, config.MAXMULTIMAP_FASTQ_PARAM)
#    aligned_bam_file = os.path.join(runconfig.output_dir, config.ALIGNED_READS_BAM_FILE)
#    aligned_log_file = os.path.join(log_dir, "bowtie_alignment.log")    
#    if all(up_to_date(aligned_bam_file, fq) for fq in runconfig.fastq_files):
#        logging.info("[SKIPPED] Alignment results exist")
#    else:    
#        logging.info("Aligning full-length reads in paired-end mode")
#        retcode = align_pe_full(runconfig.fastq_files, 
#                                bowtie_index,
#                                aligned_bam_file, 
#                                unaligned_fastq_param,
#                                maxmultimap_fastq_param,
#                                min_fragment_length=min_fragment_length,
#                                max_fragment_length=runconfig.max_fragment_length,
#                                trim5=runconfig.trim5,
#                                trim3=runconfig.trim3,
#                                library_type=runconfig.library_type,
#                                num_processors=runconfig.num_processors,
#                                fastq_format=runconfig.fastq_format,
#                                multihits=runconfig.multihits,
#                                mismatches=runconfig.mismatches,
#                                bowtie_bin=runconfig.bowtie_bin,
#                                bowtie_mode=bowtie_mode,
#                                log_file=aligned_log_file)
#        if retcode != 0:
#            logging.error("Bowtie failed with error code %d" % (retcode))    
#            sys.exit(retcode)
#    #
#    # Get insert size distribution
#    #
#    isize_dist_file = os.path.join(runconfig.output_dir, config.ISIZE_DIST_FILE)
#    isize_dist = InsertSizeDistribution()
#    if up_to_date(isize_dist_file, aligned_bam_file):
#        logging.info("[SKIPPED] Profiling insert size distribution")
#        isize_dist.from_file(open(isize_dist_file, "r"))
#    else:
#        logging.info("Profiling insert size distribution")
#        max_isize_samples = config.ISIZE_MAX_SAMPLES
#        bamfh = pysam.Samfile(aligned_bam_file, "rb")
#        isize_dist.from_bam(bamfh, min_isize=min_fragment_length, 
#                            max_isize=runconfig.max_fragment_length, 
#                            max_samples=max_isize_samples)
#        isize_dist.to_file(open(isize_dist_file, "w"))
#        bamfh.close()
#    logging.info("Insert size samples=%d mean=%f std=%f median=%d mode=%d" % 
#                 (isize_dist.n, isize_dist.mean(), isize_dist.std(), 
#                  isize_dist.percentile(50.0), isize_dist.mode()))
#    #
#    # Discordant reads alignment step
#    #
#    discordant_bam_file = os.path.join(tmp_dir, config.DISCORDANT_BAM_FILE)
#    discordant_log_file = os.path.join(log_dir, "bowtie_segmented_alignment.log")    
#    unaligned_fastq_files = [os.path.join(tmp_dir, fq) for fq in config.UNALIGNED_FASTQ_FILES]
#    # get the segments used in discordant alignment to know the effective
#    # read length used to align.  we used this to set the 'padding' during
#    # spanning read discovery
#    segments = determine_read_segments(original_read_length, 
#                                       segment_length=runconfig.segment_length, 
#                                       segment_trim=True, 
#                                       trim5=runconfig.trim5,
#                                       trim3=runconfig.trim3)
#    segmented_read_length = segments[-1][1]
#    logging.debug("Segmented alignment will use effective read length of %d" % 
#                  (segmented_read_length))
#    if all(up_to_date(discordant_bam_file, fq) for fq in runconfig.fastq_files):
#        logging.info("[SKIPPED] Discordant alignment results exist")
#    else:
#        logging.info("Aligning initially unmapped reads in single read mode")
#        align(unaligned_fastq_files, runconfig.fastq_format, bowtie_index,
#              discordant_bam_file, 
#              bowtie_bin=runconfig.bowtie_bin,
#              num_processors=runconfig.num_processors, 
#              segment_length=runconfig.segment_length,
#              segment_trim=True,
#              trim5=runconfig.trim5, 
#              trim3=runconfig.trim3, 
#              multihits=runconfig.multihits,
#              mismatches=runconfig.mismatches, 
#              bowtie_mode=bowtie_mode,
#              best_strata=runconfig.best_strata,
#              log_file=discordant_log_file)
#    #
#    # Merge paired-end reads step
#    #
#    paired_bam_file = os.path.join(tmp_dir, config.DISCORDANT_PAIRED_BAM_FILE)
#    if up_to_date(paired_bam_file, discordant_bam_file):
#        logging.info("[SKIPPED] Read pairing results exist")
#    else:
#        logging.info("Pairing aligned reads")
#        bamfh = pysam.Samfile(discordant_bam_file, "rb")
#        paired_bamfh = pysam.Samfile(paired_bam_file, "wb", template=bamfh)
#        merge_read_pairs(bamfh, paired_bamfh, 
#                         runconfig.min_fragment_length,
#                         runconfig.max_fragment_length,
#                         library_type)
#        paired_bamfh.close() 
#        bamfh.close()
#    #
#    # Find discordant reads step
#    #
#    discordant_gene_bedpe_file = \
#        os.path.join(tmp_dir, config.DISCORDANT_GENE_BEDPE_FILE)
#    discordant_genome_bedpe_file = \
#        os.path.join(tmp_dir, config.DISCORDANT_GENOME_BEDPE_FILE)
#    padding = original_read_length - segmented_read_length
#    if (up_to_date(discordant_gene_bedpe_file, paired_bam_file) and
#        up_to_date(discordant_genome_bedpe_file, paired_bam_file)):
#        logging.info("[SKIPPED] Finding discordant reads")
#    else:
#        logging.info("Finding discordant reads")
#        bamfh = pysam.Samfile(paired_bam_file, "rb")
#        find_discordant_reads(bamfh, 
#                              discordant_gene_bedpe_file,
#                              discordant_genome_bedpe_file, 
#                              gene_feature_file,
#                              max_indel_size=runconfig.max_indel_size,
#                              max_isize=runconfig.max_fragment_length,
#                              max_multihits=runconfig.multihits,
#                              library_type=library_type,
#                              padding=padding)
#        bamfh.close()
#    #
#    # Extract full sequences of the discordant reads
#    #
#    extended_discordant_gene_bedpe_file = \
#        os.path.join(tmp_dir, 
#                     config.EXTENDED_DISCORDANT_GENE_BEDPE_FILE)
#    if up_to_date(extended_discordant_gene_bedpe_file, discordant_gene_bedpe_file):
#        logging.info("[SKIPPED] Retrieving full length sequences for realignment")
#    else:
#        logging.info("Retrieving full length sequences for realignment")
#        extend_sequences(unaligned_fastq_files, 
#                         discordant_gene_bedpe_file,
#                         extended_discordant_gene_bedpe_file)
#    #
#    # Sort discordant reads
#    #
#    sorted_discordant_gene_bedpe_file = os.path.join(tmp_dir, config.SORTED_DISCORDANT_GENE_BEDPE_FILE)
#    if (up_to_date(sorted_discordant_gene_bedpe_file, extended_discordant_gene_bedpe_file)):
#        logging.info("[SKIPPED] Sorting discordant BEDPE file")
#    else:        
#        logging.info("Sorting discordant BEDPE file")
#        sort_discordant_reads(extended_discordant_gene_bedpe_file, sorted_discordant_gene_bedpe_file)        
#    #
#    # Nominate chimeras step
#    #
#    encompassing_bedpe_file = os.path.join(tmp_dir, config.ENCOMPASSING_CHIMERA_BEDPE_FILE)        
#    if (up_to_date(encompassing_bedpe_file, sorted_discordant_gene_bedpe_file)):
#        logging.info("[SKIPPED] Nominating chimeras from discordant reads")
#    else:        
#        logging.info("Nominating chimeras from discordant reads")
#        nominate_chimeras(open(sorted_discordant_gene_bedpe_file, "r"),
#                          open(encompassing_bedpe_file, "w"),
#                          gene_feature_file,                          
#                          trim=config.EXON_JUNCTION_TRIM_BP)
#    #
#    # Filter encompassing chimeras step
#    #
#    filtered_encomp_bedpe_file = \
#        os.path.join(tmp_dir,
#                     config.FILTERED_ENCOMPASSING_CHIMERA_BEDPE_FILE)
#    if (up_to_date(filtered_encomp_bedpe_file, encompassing_bedpe_file)):
#        logging.info("[SKIPPED] Filtering encompassing chimeras")
#    else:
#        logging.info("Filtering encompassing chimeras")
#        # max_isize = isize_mean + runconfig.filter_isize_stdevs*isize_std
#        filter_encompassing_chimeras(encompassing_bedpe_file,
#                                     filtered_encomp_bedpe_file,
#                                     gene_feature_file,
#                                     max_multimap=runconfig.filter_max_multimaps,
#                                     multimap_cov_ratio=runconfig.filter_multimap_ratio,
#                                     max_isize=-1,
#                                     strand_pval=runconfig.filter_strand_pval)
#    #
#    # Nominate spanning reads step
#    #
#    spanning_fastq_file = os.path.join(runconfig.output_dir, 
#                                       config.SPANNING_FASTQ_FILE)
#    if all(up_to_date(spanning_fastq_file, f) for f in unaligned_fastq_files):
#        logging.info("[SKIPPED] Preparing junction spanning reads")
#    else:
#        logging.info("Preparing junction spanning reads")
#        outfh = open(spanning_fastq_file, "w")
#        for f in unaligned_fastq_files:
#            shutil.copyfileobj(open(f), outfh)
#        outfh.close()        
#    # TODO: skip this step for now, and simply realign all the reads
##    spanning_fastq_file = os.path.join(runconfig.output_dir, config.SPANNING_FASTQ_FILE)
##    if (up_to_date(spanning_fastq_file, extended_discordant_bedpe_file) and 
##        up_to_date(spanning_fastq_file, filtered_encomp_bedpe_file)):
##        logging.info("[SKIPPED] Nominating junction spanning reads")
##    else:
##        logging.info("Nominating junction spanning reads")
##        nominate_spanning_reads(open(extended_discordant_bedpe_file, 'r'),
##                                open(filtered_encomp_bedpe_file, 'r'),
##                                open(spanning_fastq_file, 'w'))    
#    #
#    # Extract junction sequences from chimeras file
#    #        
#    ref_fasta_file = os.path.join(runconfig.index_dir, config.ALIGN_INDEX + ".fa")
#    junc_fasta_file = os.path.join(tmp_dir, config.JUNC_REF_FASTA_FILE)
#    junc_map_file = os.path.join(tmp_dir, config.JUNC_REF_MAP_FILE)
#    spanning_read_length = get_read_length(spanning_fastq_file)    
#    if (up_to_date(junc_fasta_file, filtered_encomp_bedpe_file) and
#        up_to_date(junc_map_file, filtered_encomp_bedpe_file)):        
#        logging.info("[SKIPPED] Extracting junction read sequences")
#    else:        
#        logging.info("Extracting junction read sequences")
#        bedpe_to_junction_fasta(filtered_encomp_bedpe_file, ref_fasta_file,                                
#                                spanning_read_length, 
#                                open(junc_fasta_file, "w"),
#                                open(junc_map_file, "w"))
#    #
#    # Build a bowtie index to align and detect spanning reads
#    #
#    bowtie_spanning_index = os.path.join(tmp_dir, config.JUNC_BOWTIE_INDEX)
#    bowtie_spanning_index_file = os.path.join(tmp_dir, config.JUNC_BOWTIE_INDEX_FILE)
#    if (up_to_date(bowtie_spanning_index_file, junc_fasta_file)):
#        logging.info("[SKIPPED] Building bowtie index for junction-spanning reads")
#    else:        
#        logging.info("Building bowtie index for junction-spanning reads")
#        args = [runconfig.bowtie_build_bin, junc_fasta_file, bowtie_spanning_index]
#        f = open(os.path.join(log_dir, "bowtie_build.log"), "w")
#        subprocess.call(args, stdout=f, stderr=f)
#        f.close()
#    #
#    # Align unmapped reads across putative junctions
#    #
#    junc_bam_file = os.path.join(tmp_dir, config.JUNC_READS_BAM_FILE)
#    junc_log_file = os.path.join(log_dir, "bowtie_spanning_alignment.log")        
#    if (up_to_date(junc_bam_file, bowtie_spanning_index_file) and
#        up_to_date(junc_bam_file, spanning_fastq_file)):
#        logging.info("[SKIPPED] Aligning junction spanning reads")
#    else:            
#        logging.info("Aligning junction spanning reads")
#        retcode = align_sr_full(spanning_fastq_file, 
#                                bowtie_spanning_index,
#                                junc_bam_file,
#                                trim5=runconfig.trim5,
#                                trim3=runconfig.trim3,                                 
#                                num_processors=runconfig.num_processors,
#                                fastq_format=runconfig.fastq_format,
#                                multihits=runconfig.multihits,
#                                mismatches=runconfig.mismatches,
#                                bowtie_bin=runconfig.bowtie_bin,
#                                bowtie_mode=bowtie_mode,
#                                log_file=junc_log_file)
#        if retcode != 0:
#            logging.error("Bowtie failed with error code %d" % (retcode))    
#            sys.exit(retcode)
#    #
#    # Merge spanning and encompassing read information
#    #
#    raw_chimera_bedpe_file = os.path.join(tmp_dir, config.RAW_CHIMERA_BEDPE_FILE)
#    if (up_to_date(raw_chimera_bedpe_file, junc_bam_file) and
#        up_to_date(raw_chimera_bedpe_file, junc_map_file)):
#        logging.info("[SKIPPED] Merging spanning and encompassing read alignments")
#    else:
#        logging.info("Merging spanning and encompassing read alignments")
#        merge_spanning_alignments(junc_bam_file, junc_map_file, 
#                                  raw_chimera_bedpe_file,
#                                  anchor_min=0, 
#                                  anchor_max=0,
#                                  anchor_mismatches=0)
#    #
#    # Choose best isoform for each junction
#    #
#    chimera_bedpe_file = os.path.join(tmp_dir, config.CHIMERA_BEDPE_FILE)
#    if (up_to_date(chimera_bedpe_file, raw_chimera_bedpe_file)):
#        logging.info("[SKIPPED] Filtering chimeras")
#    else:
#        logging.info("Filtering chimeras")
#        # get insert size at prob    
#        max_isize = isize_dist.percentile(runconfig.filter_isize_percentile)
#        filter_spanning_chimeras(raw_chimera_bedpe_file, 
#                                 chimera_bedpe_file,
#                                 gene_feature_file,
#                                 mate_pval=runconfig.filter_strand_pval,
#                                 max_isize=max_isize)
#    #
#    # Rank chimeras
#    #
#    ranked_chimera_bedpe_file = os.path.join(runconfig.output_dir, 
#                                             config.RANKED_CHIMERA_BEDPE_FILE)
#    if (up_to_date(ranked_chimera_bedpe_file, chimera_bedpe_file)):
#        logging.info("[SKIPPED] Ranking chimeras")
#    else:
#        logging.info("Ranking chimeras")
#        rank_chimeras(chimera_bedpe_file, ranked_chimera_bedpe_file,
#                      empirical_prob=runconfig.empirical_prob)
    #
    # Cleanup
    # 
    #shutil.rmtree(tmp_dir)
    #
    # Done
    #    
    logging.info("Finished run. Chimeras written to file %s" %
                 (ranked_chimera_bedpe_file))
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
