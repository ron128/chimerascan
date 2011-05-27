'''
Created on May 26, 2011

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

from chimerascan import pysam
from chimerascan.lib import config
from chimerascan.lib.base import FR_UNSTRANDED, FR_FIRSTSTRAND, \
    FR_SECONDSTRAND, NO_STRAND, POS_STRAND, NEG_STRAND
from chimerascan.lib.sam import get_strand, get_genomic_intervals, parse_pe_reads
from chimerascan.lib.transcriptome import build_exon_interval_trees


def get_read_strand_info(r, library_type):
    strand = get_strand(r)
    if strand == NO_STRAND:
        # try to infer strand from overlapping transcripts
        intervals = get_genomic_intervals(r)
        
    else:
        # genomic strand is known either through splice junctions
        # or because this is a strand-specific protocol
        if library_type == FR_UNSTRANDED:
            align_strand = NEG_STRAND if r.is_reverse else POS_STRAND
            is_sense = (align_strand == strand)
        elif library_type == FR_FIRSTSTRAND:
            # dUTP method: read1 is antisense and read2 is sense
            is_sense = r.is_read2
        elif library_type == FR_SECONDSTRAND:
            is_sense = r.is_read1
    return strand, is_sense
        
            


def localize_read_breakpoint(r, library_type):
    # get genomic strand of read, if available
    strand = get_strand(r)
    if strand == NO_STRAND:
        # try to infer strand from overlapping transcripts
        intervals = get_genomic_intervals(r)
    else:
        
        pass

def localize_breakpoints(bamfile, gene_file, inner_dist, library_type):
    # read transcript information and build interval trees
    exon_trees = build_exon_interval_trees(gene_file)
    # parse discordant reads
    bamfh = pysam.Samfile(bamfile, "rb")
    for pe_reads in parse_pe_reads(bamfh):
        
        pass
    bamfh.close()



def nominate_chimeras(bamfile, gene_file, max_fragment_length):
    bamfh = pysam.Samfile(bamfile, "rb")    
    bamfh.close()


def main():
    from optparse import OptionParser
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = OptionParser("usage: %prog [options] <index> <segment_length> <discordant_reads.bam> <chimeras.bedpe>")
    parser.add_option('--max-fragment-length', dest="max_fragment_length", 
                      type="int", default=1000)
    parser.add_option('--library-type', dest="library_type")
    options, args = parser.parse_args()
    index_dir = args[0]
    segment_length = int(args[1])
    bamfile = args[2]

    gene_file = os.path.join(index_dir, config.GENE_FEATURE_FILE)

    input_file = args[0]
    output_file = args[1]
    gene_file = os.path.join(options.index_dir, config.GENE_FEATURE_FILE)
    nominate_chimeras(open(input_file), 
                      open(output_file,'w'),
                      gene_file)


if __name__ == '__main__':
    main()
