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

from chimerascan.lib.seq import get_qual_conversion_func
from chimerascan.lib.base import parse_lines
import chimerascan.lib.config as config

def inspect_reads(fastq_files, output_prefix, quals):
    fqiters = [parse_lines(open(f), numlines=4) for f in fastq_files]
    outfhs = [open(output_prefix + "_%d.fq" % (x+1), "w") 
              for x in xrange(len(fastq_files))]
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
