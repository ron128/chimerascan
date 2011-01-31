'''
Created on Jan 5, 2011

@author: mkiyer
'''
'''
Created on Nov 16, 2010

@author: mkiyer
'''
import logging
import os
import xml.etree.cElementTree as etree

JOB_SUCCESS = 0
JOB_ERROR = 1

BASE_PROCESSORS = 2
MIN_SEGMENT_LENGTH = 20
ALIGN_INDEX = 'align_index'
BOWTIE_INDEX_FILE = 'align_index.1.ebwt'
GENE_REF_PREFIX = 'gene_'
GENE_FEATURE_FILE = "gene_features.txt"

ALIGNED_READS_BAM_FILE = "aligned_reads.bam"
UNALIGNED_FASTQ_PARAM = "unaligned.fq"
UNALIGNED_FASTQ_FILES = ("unaligned_1.fq", "unaligned_2.fq")
MAXMULTIMAP_FASTQ_PARAM = "maxmulti.fq"
MAXMULTIMAP_FASTQ_FILES = ("maxmulti_1.fq", "maxmulti_2.fq")
ISIZE_STATS_FILE = "isize_stats.txt"
DISCORDANT_BAM_FILE = "discordant_reads.bam"
DISCORDANT_PAIRED_BAM_FILE = "discordant_reads_paired.bam"
DISCORDANT_BEDPE_FILE = "discordant_reads.bedpe"
EXTENDED_DISCORDANT_BEDPE_FILE = "discordant_reads.extended.bedpe"
SORTED_DISCORDANT_BEDPE_FILE = "discordant_reads.srt.bedpe"
SPANNING_FASTQ_FILE = "putative_spanning_reads.fq"
ENCOMPASSING_CHIMERA_BEDPE_FILE = "encompassing_chimeras.bedpe"
EXON_JUNCTION_TRIM_BP = 10
JUNC_REF_FASTA_FILE = "spanning_juncs.fa"
JUNC_REF_MAP_FILE = "spanning_juncs.txt"
JUNC_BOWTIE_INDEX = "spanning_juncs"
JUNC_READS_BAM_FILE = "aligned_junc_reads.bam"
RAW_CHIMERA_BEDPE_FILE = "chimeras.raw.bedpe"
CHIMERA_BEDPE_FILE = "chimeras.bedpe"


translate_quals = {'solexa': 'solexa-quals',
                   'illumina': 'solexa1.3-quals',
                   'sanger': 'phred33-quals'}

class JobError(Exception):
    def __init__(self, msg):
        self.msg = msg
    def __repr__(self):
        return "'<%s(msg='%s'" % (self.__class__.__name__, self.msg)

def parse_bool(s):
    if s is None:
        return False
    return True if s.lower() in set(["t", "true", "yes", "y"]) else False

def file_exists_and_is_newer(outfile, mtime):    
    if os.path.exists(outfile):
        return os.path.getmtime(outfile) >= mtime
    return False

def file_exists_and_size_gt(outfile, filesize):    
    if os.path.exists(outfile):
        return os.path.getsize(outfile) >= filesize
    return False

def file_newer(outfile, infile):
    return file_exists_and_is_newer(outfile, os.path.getmtime(infile))

def file_newer_and_larger(outfile, infile):
    a = file_exists_and_is_newer(outfile, os.path.getmtime(infile))
    b = file_exists_and_size_gt(outfile, os.path.getsize(infile))
    return a and b

class PipelineConfig(object):
    @staticmethod
    def from_xml(xmlfile):
        config = PipelineConfig()
        tree = etree.parse(xmlfile)        
        root = tree.getroot()
        # directories
        config.output_dir = root.findtext('output_dir')
        config.tmp_dir = root.findtext('tmp_dir') 
        # fastqc tools
        config.fastqc_bin = root.findtext('fastqc_bin')
        # samtools
        config.samtools_bin = root.findtext('samtools_bin')
        # bedtools
        config.bedtools_path = root.findtext('bedtools_path')
        # bowtie
        bowtie_elem = root.find("bowtie")
        config.bowtie_path = bowtie_elem.get("path")
        config.bowtie_bin = os.path.join(config.bowtie_path, "bowtie")
        config.bowtie_build = os.path.join(config.bowtie_path, "bowtie-build")
        config.bowtie_index = bowtie_elem.findtext("index")
        config.bowtie_threads = int(bowtie_elem.findtext("processes"))
        config.multihits = int(bowtie_elem.findtext("multihits"))
        config.mismatches = int(bowtie_elem.findtext("mismatches"))
        config.seed_length = int(bowtie_elem.findtext("seed_length"))
        # genome reference
        config.ref_fasta_file = root.findtext("ref_fasta_file")
        # gene alignments 
        config.gene_bed_file = root.findtext("gene_bed_file")
        config.gene_name_file = root.findtext("gene_name_file")        
        config.gene_fasta_prefix = root.findtext("gene_fasta_prefix")
        config.insert_size_max = int(root.findtext("insert_size_max"))
        # chimera detection thresholds
        config.anchor_min = root.findtext("anchor_min")
        config.anchor_max = root.findtext("anchor_max")
        config.anchor_mismatches = root.findtext("anchor_mismatches")        
        return config

class JobConfig(object):
    @staticmethod
    def from_xml(job_file, root_output_dir):
        if not os.path.exists(root_output_dir):
            logging.error("Cannot create JobConfig because output_dir %s does not exist" % (root_output_dir))
        j = JobConfig()
        tree = etree.parse(job_file)    
        root = tree.getroot()
        output_dir_elem = root.find("output_dir")
        j.dst_output_dir = output_dir_elem.text
        j.remote = parse_bool(output_dir_elem.get("remote"))
        if j.remote:
            j.remote_ip = output_dir_elem.get("ip")
        else:
            j.remote_ip = None
        j.name = root.get("name")
        j.output_dir = os.path.join(root_output_dir, j.name)
        j.species = root.findtext("species")
        j.fragment_layout = root.findtext("fragment_layout")
        j.fragment_length_mean = int(root.findtext("fragment_length_mean"))
        j.fastq_format = translate_quals[root.findtext("quality_scores")]    
        j.read_length = int(root.findtext("read_length"))
        src_fastq_files = {}
        dst_fastq_files = {}
        fastq_files = {}
        for mate_elem in root.findall("fastq_files/*"):
            mate = int(mate_elem.get("mate"))
            src_fastq_files[mate] = mate_elem.text
            ext = os.path.splitext(mate_elem.text)[1]
            dst_fastq_files[mate] = os.path.join(j.output_dir, "mate%d%s" % (mate,ext))
            fastq_files[mate] = os.path.join(j.output_dir, "mate%d.fq" % (mate))
        j.src_fastq_files = [src_fastq_files[mate] for mate in xrange(len(src_fastq_files))]
        j.dst_fastq_files = [dst_fastq_files[mate] for mate in xrange(len(dst_fastq_files))]
        j.fastq_files = [fastq_files[mate] for mate in xrange(len(fastq_files))]
        # alignment step
        j.chimerascan_dir = os.path.join(j.output_dir, "chimerascan")
        j.discordant_bam_file = os.path.join(j.chimerascan_dir, DISCORDANT_READS_FILE)
        j.expression_file = os.path.join(j.chimerascan_dir, EXPRESSION_FILE)
        # chimera nomination
        j.chimera_bedpe_file = os.path.join(j.chimerascan_dir, "filtered_chimeras.bedpe.txt")
        # chimera fasta file
        j.chimera_fasta_file = os.path.join(j.chimerascan_dir, "chimericjuncs.fa")
        j.chimera_mapping_file = os.path.join(j.chimerascan_dir, "chimericjunc_mapping.txt")
        # spanning bowtie index
        j.bowtie_chimera_index = os.path.splitext(j.chimera_fasta_file)[0]
        # spanning read output files
        j.spanning_bowtie_output_files = ["mate%d_spanning_bowtie.txt" % (mate) for mate in xrange(len(fastq_files))]
        # chimera output file
        j.spanning_chimera_file = os.path.join(j.chimerascan_dir, "chimeras.bedpe.txt")
        return j
