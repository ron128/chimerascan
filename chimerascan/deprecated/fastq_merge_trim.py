'''
Created on May 23, 2011

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
import argparse

from chimerascan.lib import config

def parse_fastq(line_iter):
    with line_iter:
        while True:
            lines = [line_iter.next().rstrip() for x in xrange(4)]
            yield lines

def trim_and_merge_fastq(infiles, outfile, segment_length):
    fqiters = [parse_fastq(open(f)) for f in infiles]    
    if outfile == "-":
        outfh = sys.stdout
    else:
        outfh = open(outfile, "w")
    try:
        while True:
            pe_lines = [fqiter.next() for fqiter in fqiters]
            for readnum,lines in enumerate(pe_lines):
                seqlen = len(lines[1])
                if seqlen > segment_length:
                    # encode a '0' or '1' as the first character of the line
                    lines[0] = "@%d%s" % (readnum, lines[0][1:])
                    lines[1] = lines[1][:segment_length]
                    lines[3] = lines[3][:segment_length]
                print >>outfh, '\n'.join(lines)
    except StopIteration:
        pass
    if outfile != "-":
        outfh.close()
    return config.JOB_SUCCESS

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--segment-length", type="int", dest="segment_length", default=config.MIN_SEGMENT_LENGTH)
    parser.add_argument("output_file")    
    parser.add_argument("fastq_files", nargs=2)
    args = parser.parse_args()
    return trim_and_merge_fastq(args.output_file, (args.read1, args.read2), args.segment_length)

if __name__ == '__main__':
    sys.exit(main())
