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

def parse_fastq(line_iter):
    with line_iter:
        while True:
            lines = [line_iter.next().rstrip() for x in xrange(4)]
            yield lines

def trim_and_merge_fastq(infiles, outfile, trim5, segment_length):
    total_length = trim5 + segment_length
    fqiters = [parse_fastq(open(f)) for f in infiles]    
    if outfile == "-":
        outfh = sys.stdout
    else:
        outfh = open(outfile, "w")
    try:
        while True:
            pe_lines = [fqiter.next() for fqiter in fqiters]
            for lines in pe_lines:
                seqlen = len(lines[1])
                if seqlen > total_length:
                    lines[1] = lines[1][trim5:total_length]
                    lines[3] = lines[3][trim5:total_length]
                print >>outfh, '\n'.join(lines)
    except StopIteration:
        pass
    if outfile != "-":
        outfh.close()

def main():
    from optparse import OptionParser
    parser = OptionParser("usage: %prog [options] <in1.fq> <in2.fq> <out.fq>")
    parser.add_option("--trim5", type="int", dest="trim5", default=0)
    parser.add_option("--segment-length", type="int", dest="segment_length", default=25)
    options, args = parser.parse_args()
    trim_and_merge_fastq(args[:2], args[2], options.trim5, options.segment_length)

if __name__ == '__main__':
    main()
