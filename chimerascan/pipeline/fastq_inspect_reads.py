'''
Created on Jul 14, 2011

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
import logging
import gzip
import bz2
import zipfile
import os

from chimerascan.lib.seq import get_qual_conversion_func
from chimerascan.lib.base import parse_lines
import chimerascan.lib.config as config

def detect_format(fastq_files):
    if all(f.endswith(".gz") for f in fastq_files):
        return "gz"
    elif all(f.endswith(".bz2") for f in fastq_files):
        return "bz2"
    elif all(f.endswith(".zip") for f in fastq_files):
        return "zip"
    else:
        return "txt"

def open_compressed(fastq_files):
    compression_format = detect_format(fastq_files)
    if compression_format == "gz":
        filehandles = [gzip.open(f, "r") for f in fastq_files]
    elif compression_format == "bz2":
        filehandles = [bz2.BZ2File(f, "r") for f in fastq_files]
    elif compression_format == "zip":
        filehandles = [zipfile.ZipFile(f, "r") for f in fastq_files]
    else:
        filehandles = [open(f, "r") for f in fastq_files]
    return filehandles

def detect_read_lengths(fastq_files):
    filehandles = open_compressed(fastq_files)
    tags = [f.next() for f in filehandles]
    seqs = [f.next() for f in filehandles]
    return [len(s) for s in seqs]

def inspect_reads(fastq_files, output_prefix, quals):
    """
    uncompresses reads, renames reads, and converts quality scores 
    to 'sanger' format
    """
    # setup file iterators
    filehandles = open_compressed(fastq_files)
    fqiters = [parse_lines(f, numlines=4) for f in filehandles]
    output_files = [(output_prefix + "_%d.fq" % (x+1)) 
                    for x in xrange(len(fastq_files))]
    outfhs = [open(f, "w") for f in output_files]
    qual_func = get_qual_conversion_func(quals)
    linenum = 0    
    try:
        while True:
            pelines = [it.next() for it in fqiters]
            for i,lines in enumerate(pelines):
                # rename read using line number
                lines[0] = "@%d/%d" % (linenum,i+1)
                # ignore redundant header
                lines[2] = "+"
                # convert quality score to sanger
                lines[3] = qual_func(lines[3])
                print >>outfhs[i], '\n'.join(lines)
            linenum += 1
    except StopIteration:
        pass
    except:
        logging.error("Unexpected error during FASTQ file processing")
        for f in output_files:
            if os.path.exists(f):
                os.remove(f)
        return config.JOB_ERROR
    for fh in filehandles:
        fh.close()
    logging.debug("Inspected %d fragments" % (linenum))
    return config.JOB_SUCCESS

def main():
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")    
    from optparse import OptionParser
    parser = OptionParser("usage: %prog [options] <outprefix> <in1.fq> <in2.fq>")
    parser.add_option("--quals", dest="quals", choices=["sanger", "solexa", "illumina"], 
                      default="sanger")
    options, args = parser.parse_args()
    if len(args) < 2:
        parser.error("must specify output prefix and at least one fastq file")
    output_prefix = args[0]
    fastq_files = args[1:]
    inspect_reads(fastq_files, output_prefix, options.quals)

if __name__ == '__main__':
    main()
