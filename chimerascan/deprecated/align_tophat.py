'''
Created on Apr 27, 2011

@author: mkiyer

chimerascan: chimeric transcript discovery using RNA-seq

Copyright (C) 2011 Matthew Iyer, Christopher Maher

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
import logging
import subprocess
import sys
import os

from chimerascan.lib.seq import parse_fastq
from chimerascan.lib.base import make_temp
from fastq_merge_trim import trim_and_merge_fastq

translate_quals = {'solexa': '--solexa-quals',
                   'illumina': '--solexa1.3-quals'}

def align_tophat(fastq_files, tophat_bin, bowtie_index, output_dir,
                 read_length, insert_size, insert_size_stdev,
                 num_processors, library_type, quals, multihits,
                 mismatches, segment_length, raw_juncs_file, 
                 log_file=None):
    # compute mate inner dist
    mate_inner_dist = insert_size - 2*read_length
    if mate_inner_dist < 0:
        logging.warning("Inner distance between reads is estimated to be "
                        "negative, indicating that the read alignments are "
                        "likely to overlap")
    # setup arguments for tophat
    args = [tophat_bin, 
            "-o", output_dir,
            "-r", mate_inner_dist,
            "--mate-std-dev", insert_size_stdev,
            "-g", multihits,
            "-p", num_processors,
            "--library-type", library_type,
            "--segment-mismatches", mismatches,
            "--segment-length", segment_length]
    if raw_juncs_file is not None:
        args.extend(["-j", raw_juncs_file])
    if quals in translate_quals:
        args.append(translate_quals[quals])
    # positional arguments
    args.append(bowtie_index)
    args.extend(fastq_files)
    # kickoff process
    args = map(str, args)
    logging.debug("Tophat alignment args: %s" % (' '.join(args)))
    if log_file is not None:
        logfh = open(log_file, "w")
    else:
        logfh = sys.stderr
    retcode = subprocess.call(args, stderr=logfh)
    if log_file is not None:
        logfh.close()
    if retcode != 0:
        logging.error("Tophat exited with error code '%d'" % (retcode))
    return retcode

def realign_tophat_sr(fastq_file, tophat_bin, bowtie_index, output_dir,
                      num_processors, library_type, quals, multihits,
                      mismatches, segment_length, raw_juncs_file, 
                      log_file=None):
    # setup arguments for tophat
    args = [tophat_bin, 
            "-o", output_dir,
            "-g", multihits,
            "-p", num_processors,
            "--library-type", library_type,
            "--segment-mismatches", mismatches,
            "--segment-length", segment_length,
            "--no-novel-juncs"]
    if raw_juncs_file is not None:
        args.extend(["-j", raw_juncs_file])
    if quals in translate_quals:
        args.append(translate_quals[quals])
    # positional arguments
    args.append(bowtie_index)
    args.append(fastq_file)
    # kickoff process
    args = map(str, args)
    logging.debug("Tophat alignment args: %s" % (' '.join(args)))
    if log_file is not None:
        logfh = open(log_file, "w")
    else:
        logfh = sys.stderr
    retcode = subprocess.call(args, stderr=logfh)
    if log_file is not None:
        logfh.close()
    if retcode != 0:
        logging.error("Tophat exited with error code '%d'" % (retcode))
    return retcode

