'''
Created on Jul 6, 2011

@author: mkiyer
'''
'''
Created on Jun 3, 2011

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

from chimerascan.lib.seq import parse_fastq_record, get_qual_conversion_func

def inspect_reads(fastq_files, output_fastq_files, qual_format):
    fqiters = [parse_fastq_record(open(f)) for f in fastq_files]
    outfhs = [open(f, "w") for f in output_fastq_files]    
    qual_func = get_qual_conversion_func(qual_format)
    try:
        while True:
            fqrecs = [it.next() for it in fqiters]
            for i,fqrec in enumerate(fqrecs):
                fqrec.qual = qual_func(fqrec.qual)
                print >>outfhs[i], fqrec.to_string()
    except StopIteration:
        pass

def main():
    from optparse import OptionParser
    parser = OptionParser("usage: %prog [options] <in1.fq> <out1.fq> [<in2.fq> <out2.fq>]")
    parser.add_option("--quals", dest="quals", choices=["sanger", "solexa", "illumina"], 
                      default="sanger")
    options, args = parser.parse_args()
    if len(args) < 2:
        parser.error("must specify at least one fastq file and one output file")
    fastq_files = [args[0]]
    output_files = [args[1]]
    if len(args) >= 4:
        fastq_files.append(args[2])
        output_files.append(args[3])
    inspect_reads(fastq_files, output_files, options.quals)

if __name__ == '__main__':
    main()
