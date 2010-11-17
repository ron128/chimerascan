'''
Created on Nov 16, 2010

@author: mkiyer
'''
import logging
import os
import xml.etree.cElementTree as etree

JOB_SUCCESS = 0
JOB_ERROR = 1

FASTQC_DATA_FILE = 'fastqc_data.txt'
FASTQC_REPORT_FILE = 'fastqc_report.html'

DISCORDANT_READS_FILE = "discordant_reads.bam"
EXPRESSION_FILE = "expression.txt"

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
        # bowtie
        bowtie_elem = root.find("bowtie")
        config.bowtie_bin = bowtie_elem.get("bin")
        config.bowtie_index = bowtie_elem.findtext("index")
        config.bowtie_threads = int(bowtie_elem.findtext("processes"))
        config.multihits = int(bowtie_elem.findtext("multihits"))
        config.mismatches = int(bowtie_elem.findtext("mismatches"))
        config.seed_length = int(bowtie_elem.findtext("seed_length"))
        # gene alignments 
        config.gene_bed_file = root.findtext("gene_bed_file")
        config.gene_fasta_prefix = root.findtext("gene_fasta_prefix")
        config.insert_size_max = int(root.findtext("insert_size_max"))
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
        # chimera scan
        j.chimerascan_dir = os.path.join(j.output_dir, "chimerascan")
        j.discordant_bam_file = os.path.join(j.chimerascan_dir, DISCORDANT_READS_FILE)
        j.expression_file = os.path.join(j.chimerascan_dir, EXPRESSION_FILE)
        return j
