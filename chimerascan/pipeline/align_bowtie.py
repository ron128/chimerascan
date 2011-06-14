'''
Created on Jun 1, 2011

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
import sys
import os
import logging
import subprocess

from chimerascan.lib.base import get_read_length, LibraryTypes
from chimerascan.lib.seq import SANGER_FORMAT, SOLEXA_FORMAT, ILLUMINA_FORMAT
from chimerascan.lib.config import MIN_SEGMENT_LENGTH, JOB_SUCCESS

translate_quals = {SOLEXA_FORMAT: 'solexa-quals',
                   ILLUMINA_FORMAT: 'solexa1.3-quals',
                   SANGER_FORMAT: 'phred33-quals'}

def translate_library_type(library_type):
    """
    returns the bowtie library type option '--fr' or '--ff' corresponding
    to the first two characters of the library type string
    """
    return library_type[0:2]

_sam2bam_script = os.path.join(os.path.dirname(__file__), "sam2bam.py")
_fastq_trim_script = os.path.join(os.path.dirname(__file__), "fastq_merge_trim.py")
_discordant_script = os.path.join(os.path.dirname(__file__), "find_discordant_reads.py")

def align_pe(fastq_files, 
             bowtie_index,
             output_bam_file, 
             unaligned_fastq_param=None,
             maxmultimap_fastq_param=None,
             min_fragment_length=0,
             max_fragment_length=1000,
             trim5=0,
             trim3=0,
             library_type=LibraryTypes.FR_UNSTRANDED,
             num_processors=1, 
             quals=SANGER_FORMAT,
             multihits=100, 
             mismatches=2, 
             bowtie_bin="bowtie", 
             bowtie_args=None,
             log_file=None,
             keep_unmapped=False):
    read_length = get_read_length(fastq_files[0])
    args = [bowtie_bin, "-q", "-S", 
            "-p", str(num_processors),
            "--%s" % translate_quals[quals],
            "-k", str(multihits),
            "-m", str(multihits),
            "-n", str(mismatches),
            "-l", str(read_length),
            "--minins", min_fragment_length,
            "--maxins", max_fragment_length,
            "--trim5", trim5,
            "--trim3", trim3,
            "--%s" % translate_library_type(library_type)]
    if unaligned_fastq_param is not None:
        args.extend(["--un", unaligned_fastq_param])
    if maxmultimap_fastq_param is not None:
        args.extend(["--max", maxmultimap_fastq_param])    
    if bowtie_args is not None:        
        args.extend(bowtie_args.split())
    args += [bowtie_index, 
             "-1", fastq_files[0],
             "-2", fastq_files[1]]
    args = map(str, args)
    logging.debug("Bowtie alignment args: %s" % (' '.join(args)))
    # setup logging
    if log_file is not None:
        logfh = open(log_file, "w")
    else:
        logfh = None
    aln_p = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=logfh)
    # pipe the bowtie SAM output to a filter that writes BAM format
    args = [sys.executable, _sam2bam_script, 
            "--multihits", str(multihits),
            "--quals", quals]
    if keep_unmapped:
        args.append("--un")
    args.extend([output_bam_file, "-"])
    args.extend(fastq_files)
    logging.debug("SAM to BAM converter args: %s" % (' '.join(args)))
    retcode = subprocess.call(args, stdin=aln_p.stdout, stderr=logfh)    
    if logfh is not None:
        logfh.close()
    if retcode != 0:
        aln_p.terminate()
        return retcode
    return aln_p.wait()

def align_sr(fastq_file, 
             bowtie_index,
             output_bam_file, 
             unaligned_fastq_param=None,
             maxmultimap_fastq_param=None,
             trim5=0,
             trim3=0,
             library_type=LibraryTypes.FR_UNSTRANDED,
             num_processors=1, 
             quals=SANGER_FORMAT,
             multihits=100, 
             mismatches=2, 
             bowtie_bin="bowtie", 
             bowtie_args=None,
             log_file=None,
             keep_unmapped=False):
    args = [bowtie_bin, "-q", "-S", 
            "-p", str(num_processors),
            "--%s" % translate_quals[quals],
            "-k", str(multihits),
            "-m", str(multihits),
            "-n", str(mismatches),
            "--trim5", trim5,
            "--trim3", trim3,
            "--%s" % translate_library_type(library_type)]
    if unaligned_fastq_param is not None:
        args.extend(["--un", unaligned_fastq_param])
    if maxmultimap_fastq_param is not None:
        args.extend(["--max", maxmultimap_fastq_param])    
    if bowtie_args is not None:        
        args.extend(bowtie_args.split())
    args += [bowtie_index, fastq_file]
    args = map(str, args)
    logging.debug("Bowtie alignment args: %s" % (' '.join(args)))
    # setup logging
    if log_file is not None:
        logfh = open(log_file, "w")
    else:
        logfh = None
    #
    # Start bowtie alignment process
    # 
    aln_p = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=logfh)
    #
    # Fix alignment ordering and convert to BAM, also extend sequences
    # back to full length by adding padding to CIGAR string
    #
    args = [sys.executable, _sam2bam_script, 
            "--multihits", str(multihits),
            "--quals", quals,
            "--pesr"]
    if keep_unmapped:
        args.append("--un")
    args.extend([output_bam_file, "-"])
    args.append(fastq_file)    
    logging.debug("SAM to BAM converter args: %s" % (' '.join(args)))
    fix_p = subprocess.Popen(args, stdin=aln_p.stdout, stderr=logfh)
    # wait for processes to complete
    retcode1 = fix_p.wait()
    retcode2 = aln_p.wait()
    # end logging
    if logfh is not None:
        logfh.close()
    return retcode1 or retcode2


def trim_align_pe_sr(fastq_files,
                     bowtie_index,
                     output_bam_file,
                     unaligned_fastq_param=None,
                     maxmultimap_fastq_param=None,
                     trim5=0,
                     library_type=LibraryTypes.FR_UNSTRANDED,
                     num_processors=1, 
                     quals=SANGER_FORMAT,
                     multihits=100, 
                     mismatches=2, 
                     bowtie_bin="bowtie", 
                     bowtie_args=None,
                     log_file=None,
                     segment_length=25,
                     keep_unmapped=False):
    # setup logging
    if log_file is not None:
        logfh = open(log_file, "w")
    else:
        logfh = None
    #
    # Merge paired-end reads into single fastq file
    #
    args = [sys.executable, _fastq_trim_script, 
            "--trim5", str(trim5), 
            "--segment-length", str(segment_length)]
    args.extend(fastq_files)
    args.append("-")
    logging.debug("FASTQ trimming args: %s" % (' '.join(args)))
    trim_p = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=logfh)
    #
    # Align the trimmed reads
    #
    args = [bowtie_bin, "-q", "-S", 
            "-p", str(num_processors),
            "--tryhard",
            "--%s" % translate_quals[quals],
            "-k", str(multihits),
            "-m", str(multihits),
            "-n", str(mismatches),
            "-l", str(segment_length),
            "--%s" % translate_library_type(library_type)]
    if unaligned_fastq_param is not None:
        args.extend(["--un", unaligned_fastq_param])
    if maxmultimap_fastq_param is not None:
        args.extend(["--max", maxmultimap_fastq_param])            
    if bowtie_args is not None:        
        args.extend(bowtie_args.split())
    args += [bowtie_index, "-"]
    logging.debug("Alignment args: %s" % (' '.join(args)))
    aln_p = subprocess.Popen(args, stdin=trim_p.stdout, 
                             stdout=subprocess.PIPE,
                             stderr=logfh)
    #
    # Fix alignment ordering and convert to BAM, also extend sequences
    # back to full length by adding padding to CIGAR string
    #
    args = [sys.executable, _sam2bam_script, 
            "--multihits", str(multihits),
            "--quals", quals,
            "--pesr", 
            "--softclip"] 
    if keep_unmapped:
        args.append("--un")
    args.extend([output_bam_file, "-"])
    args.extend(fastq_files)
    logging.debug("SAM to BAM converter args: %s" % (' '.join(args)))
    fix_p = subprocess.Popen(args, stdin=aln_p.stdout, stderr=logfh)
    # wait for processes to complete
    fix_p.wait()
    aln_p.wait()
    trim_p.wait()
    # end logging
    if logfh is not None:
        logfh.close()

