'''
Created on Jul 1, 2011

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
import sys
import operator

from chimerascan.lib.chimera import Chimera, CHIMERA_SEP
from chimerascan.lib import config
from chimerascan.lib.base import make_temp
from chimerascan.lib.gene_to_genome import build_gene_to_genome_map, gene_to_genome_pos

def write_output(input_file, output_file, index_dir):
    gene_file = os.path.join(index_dir, config.GENE_FEATURE_FILE)
    # build a lookup table to get genome coordinates from transcript 
    # coordinates
    tx_genome_map = build_gene_to_genome_map(open(gene_file))
    # read chimera data into memory
    # TODO: could explode if many chimeras exist
    lines = []
    for c in Chimera.parse(open(input_file)):
        chrom5p,strand5p,pos5p = gene_to_genome_pos(c.partner5p.tx_name, c.partner5p.end-1, tx_genome_map)
        chrom3p,strand3p,pos3p = gene_to_genome_pos(c.partner3p.tx_name, c.partner3p.start, tx_genome_map)
        fields = [c.partner5p.tx_name, c.partner5p.start, c.partner5p.end,
                  c.partner3p.tx_name, c.partner3p.start, c.partner3p.end,
                  CHIMERA_SEP.join([c.partner5p.gene_name, 
                                    c.partner3p.gene_name]),
                  c.get_weighted_cov(),
                  c.partner5p.strand, c.partner3p.strand,
                  c.name,
                  "%s:%d" % (chrom5p, pos5p),
                  "%s:%d" % (chrom3p, pos3p),
                  c.chimera_type, c.distance, 
                  c.get_total_unique_reads(), 
                  c.get_unique_spanning_reads()]
        lines.append(fields)
    # sort
    lines = sorted(lines, key=operator.itemgetter(16, 7, 15), reverse=True)    
    f = open(output_file, "w")
    print >>f, '\t'.join(['#gene5p', 'start5p', 'end5p', 'gene3p', 
                          'start3p', 'end3p', 'name', 'multimap_weighted_cov', 
                          'strand5p', 'strand3p', 'chimera_id',
                          'breakpoint5p', 'breakpoint3p',
                          'type', 'distance', 
                          'total_unique_frags', 
                          'spanning_unique_frags'])
    for fields in lines:
        print >>f, '\t'.join(map(str, fields))
    f.close()
    return config.JOB_SUCCESS

def main():
    from optparse import OptionParser
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = OptionParser("usage: %prog [options] <index_dir> <in.txt> <out.txt>")
    options, args = parser.parse_args()
    index_dir = args[0]
    input_file = args[1]
    output_file = args[2]
    return write_output(input_file, output_file, index_dir)

if __name__ == "__main__":
    sys.exit(main())
