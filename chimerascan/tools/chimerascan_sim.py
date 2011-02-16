#!/usr/bin/env python
'''
Created on Jan 8, 2011

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

from lib.seq import DNA_reverse_complement
from lib.base import parse_library_type

def parse_fasta(line_iter):    
    line = line_iter.next()
    assert line.startswith('>')
    tag = line.rstrip()[1:]
    seq = ''
    for line in line_iter:
        if line.startswith('>'):
            yield tag, seq
            tag = line.rstrip()[1:]
            seq = ''
        else:
            seq += line.rstrip()
    yield tag, seq

def generate_library_in_silico(num_reads,
                               reference_fasta_file,
                               isize_mean=200,
                               isize_std=20,
                               frac_genes=0.9,
                               frac_genome=0.1):
    '''
    Usage:   wgsim [options] <in.ref.fa> <out.read1.fq> <out.read2.fq>
    
    Options: -e FLOAT      base error rate [0.020]
             -d INT        outer distance between the two ends [500]
             -s INT        standard deviation [50]
             -N INT        number of read pairs [1000000]
             -1 INT        length of the first read [70]
             -2 INT        length of the second read [70]
             -r FLOAT      rate of mutations [0.0010]
             -R FLOAT      fraction of indels [0.10]
             -X FLOAT      probability an indel is extended [0.30]
             -c            generate reads in color space (SOLiD reads)
             -C            show mismatch info in comment rather than read name
             -h            haplotype mode
    '''
    
    from optparse import OptionParser
    
    parser = OptionParser("usage: %prog [options] <seq.fa> <outprefix>")
    parser.add_option("--rlen", dest="read_length", type="int", default=50)
    parser.add_option("--isize", dest="insert_size", type="int", default=200)
    parser.add_option("--library", dest="library_type", default="fr")
    options, args = parser.parse_args()
    # extract command line arguments
    fasta_file = args[0]
    outprefix = args[1]
    library_type = parse_library_type(options.library_type)
    
    fasta_records = list(parse_fasta(open(fasta_file)))
    tag, seq = fasta_records[0]
    tag = tag.split()[0]
    read_length = options.read_length
    insert_size = options.insert_size
    assert insert_size > read_length
    read1 = (0, read_length)
    read2 = (insert_size - read_length, insert_size)

    fh1 = open(outprefix + '1.fq', 'w')
    fh2 = open(outprefix + '2.fq', 'w')
    for end in xrange(insert_size, len(seq)):
        start = end - insert_size        
        start1,end1 = start + read1[0], start + read1[1]
        start2,end2 = start + read2[0], start + read2[1]        
        seq1 = seq[start1:end1]
        if library_type[0] == 1:
            seq1 = DNA_reverse_complement(seq1)
        seq2 = seq[start2:end2]
        if library_type[1] == 1:
            seq2 = DNA_reverse_complement(seq2)        
        thistag = '%s:%d-%d' % (tag, start1, end2)
        print >>fh1, '@%s/1\n%s\n+%s/1\n%s' % (thistag, seq1, thistag, 'I'*read_length)        
        print >>fh2, '@%s/2\n%s\n+%s/2\n%s' % (thistag, seq2, thistag, 'I'*read_length)
    fh1.close()
    fh2.close()


def main():
    from optparse import OptionParser
    parser = OptionParser("usage: %prog [options] <seq.fa> <outprefix>")
    parser.add_option("--rlen", dest="read_length", type="int", default=50)
    parser.add_option("--isize", dest="insert_size", type="int", default=200)
    parser.add_option("--library", dest="library_type", default="fr")
    options, args = parser.parse_args()
    # extract command line arguments
    fasta_file = args[0]
    outprefix = args[1]
    library_type = parse_library_type(options.library_type)
    
    fasta_records = list(parse_fasta(open(fasta_file)))
    tag, seq = fasta_records[0]
    tag = tag.split()[0]
    read_length = options.read_length
    insert_size = options.insert_size
    assert insert_size > read_length
    read1 = (0, read_length)
    read2 = (insert_size - read_length, insert_size)

    fh1 = open(outprefix + '1.fq', 'w')
    fh2 = open(outprefix + '2.fq', 'w')
    for end in xrange(insert_size, len(seq)):
        start = end - insert_size        
        start1,end1 = start + read1[0], start + read1[1]
        start2,end2 = start + read2[0], start + read2[1]        
        seq1 = seq[start1:end1]
        if library_type[0] == 1:
            seq1 = DNA_reverse_complement(seq1)
        seq2 = seq[start2:end2]
        if library_type[1] == 1:
            seq2 = DNA_reverse_complement(seq2)        
        thistag = '%s:%d-%d' % (tag, start1, end2)
        print >>fh1, '@%s/1\n%s\n+%s/1\n%s' % (thistag, seq1, thistag, 'I'*read_length)        
        print >>fh2, '@%s/2\n%s\n+%s/2\n%s' % (thistag, seq2, thistag, 'I'*read_length)
    fh1.close()
    fh2.close()

#@PATHBIO-SOLEXA2_30TUEAAXX:1:3:1574:973/1
#ATAAGTACAGTCTATTGATCCGTGATGCATGCTGTACATGAATTTGGTGNNNN
#+PATHBIO-SOLEXA2_30TUEAAXX:1:3:1574:973/1
#bbaabbaaaabbaaaaaaaWaaaaaaZaaaaaaaaaaaaa_[_a_X_^UEEEE
                                          
if __name__ == '__main__':
    main()
