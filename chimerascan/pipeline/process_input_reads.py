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
import os
import argparse

from chimerascan.lib.seq import get_qual_conversion_func
from chimerascan.lib.base import parse_lines, open_compressed
import chimerascan.lib.config as config

def process_input_reads(fastq_files, output_prefix, quals, trim5, trim3):
    """
    uncompresses reads, renames reads, and converts quality scores 
    to 'sanger' format
    """
    # setup file iterators for input fastq files
    infhs = [open_compressed(f) for f in fastq_files]
    fqiters = [parse_lines(f, numlines=4) for f in infhs]
    # setup output files
    output_files = [(output_prefix + "_%d.fq" % (x+1)) 
                    for x in xrange(len(fastq_files))]
    outfhs = [open(f, "w") for f in output_files]
    read_name_file = output_prefix + ".txt"
    read_name_fh = open(read_name_file, 'w')
    # get quality score conversion function
    qual_func = get_qual_conversion_func(quals)
    linenum = 1
    try:
        while True:
            pelines = [it.next() for it in fqiters]
            # get read1 first line of fq record, and remove "@" symbol
            read1_name = pelines[0][0][1:]
            # remove whitespace and/or read number tags /1 or /2
            read1_name = read1_name.split()[0].split("/")[0]
            # write to read name database
            print >>read_name_fh, read1_name
            # convert reads
            for i,lines in enumerate(pelines):
                # rename read using line number
                lines[0] = "@%d/%d" % (linenum,i+1)
                # ignore redundant header
                lines[2] = "+"
                # trim read
                total_length = len(lines[1])
                pos3p = max(trim5+1, total_length - trim3)
                lines[1] = lines[1][trim5:pos3p]
                lines[3] = lines[3][trim5:pos3p]
                # convert quality score to sanger
                lines[3] = qual_func(lines[3])
                print >>outfhs[i], '\n'.join(lines)
            linenum += 1
    except StopIteration:
        pass
    except:
        logging.error("Unexpected error during FASTQ file processing")
        for fh in outfhs:
            fh.close()
        read_name_fh.close()
        for f in output_files:
            if os.path.exists(f):
                os.remove(f)
        if os.path.exists(read_name_file):
            os.remove(read_name_file)
        return config.JOB_ERROR
    # cleanup
    for fh in infhs:
        fh.close()
    for fh in outfhs:
        fh.close()
    read_name_fh.close()
    logging.debug("Inspected %d fragments" % (linenum))
    return config.JOB_SUCCESS

def main():
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")    
    parser = argparse.ArgumentParser()
    parser.add_argument("--quals", dest="quals", choices=["sanger", "solexa", "illumina"], 
                        default="sanger")
    parser.add_argument("--trim5", dest="trim5", type=int, default=0)
    parser.add_argument("--trim3", dest="trim3", type=int, default=0)
    parser.add_argument("output_prefix")
    parser.add_argument("fastq_files", nargs="+")    
    args = parser.parse_args()
    process_input_reads(args.fastq_files, args.output_prefix, args.quals, args.trim5, args.trim3)

if __name__ == '__main__':
    main()
