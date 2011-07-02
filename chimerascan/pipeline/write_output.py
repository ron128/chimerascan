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
import collections

from chimerascan.lib.chimera import Chimera, CHIMERA_SEP
from chimerascan.lib import config
from chimerascan.lib.base import make_temp
from chimerascan.lib.gene_to_genome import build_gene_to_genome_map, \
    build_tx_cluster_map, gene_to_genome_pos

def get_chimera_groups(input_file, gene_file):
    # build a lookup table to get gene clusters from transcript name    
    tx_cluster_map = build_tx_cluster_map(open(gene_file))
    # build a lookup table to get genome coordinates from transcript 
    # coordinates
    tx_genome_map = build_gene_to_genome_map(open(gene_file))
    # group chimeras in the same genomic cluster with the same
    # breakpoint
    cluster_chimera_dict = collections.defaultdict(lambda: [])
    for c in Chimera.parse(open(input_file)):
        # get cluster of overlapping genes
        cluster5p = tx_cluster_map[c.partner5p.tx_name]
        cluster3p = tx_cluster_map[c.partner3p.tx_name]
        # get genomic positions of breakpoints
        coord5p = gene_to_genome_pos(c.partner5p.tx_name, c.partner5p.end-1, tx_genome_map)
        coord3p = gene_to_genome_pos(c.partner3p.tx_name, c.partner3p.start, tx_genome_map)
        # add to dictionary
        cluster_chimera_dict[(cluster5p,cluster3p,coord5p,coord3p)].append(c)
    for key,chimeras in cluster_chimera_dict.iteritems():
        yield key,chimeras

def write_output(input_file, output_file, index_dir):
    gene_file = os.path.join(index_dir, config.GENE_FEATURE_FILE)
    # group chimera isoforms together
    # TODO: requires reading all chimeras into memory
    lines = []
    for key,chimeras in get_chimera_groups(input_file, gene_file):
        cluster5p,cluster3p,coord5p,coord3p = key
        chrom5p,strand5p,pos5p = coord5p
        chrom3p,strand3p,pos3p = coord3p
        txs5p = ",".join(set(c.partner5p.tx_name for c in chimeras))
        txs3p = ",".join(set(c.partner3p.tx_name for c in chimeras))
        genes5p = ",".join(set(c.partner5p.gene_name for c in chimeras))
        genes3p = ",".join(set(c.partner3p.gene_name for c in chimeras))
        c = chimeras[0]
        fields = [chrom5p, pos5p, strand5p, 
                  chrom3p, pos3p, strand3p,
                  txs5p, txs3p, genes5p, genes3p,
                  c.chimera_type, c.distance,
                  c.get_weighted_cov(),
                  c.get_total_unique_reads(),
                  c.num_encomp_frags,
                  c.get_unique_spanning_reads()]
        lines.append(fields)
    # sort
    lines = sorted(lines, key=operator.itemgetter(15, 14, 12), reverse=True)    
    f = open(output_file, "w")
    print >>f, '\t'.join(['#chrom5p', 'breakpoint_pos_5p', 'strand5p',
                          'chrom3p', 'breakpoint_pos_3p', 'strand3p',
                          'transcript_ids_5p', 'transcript_ids_3p',
                          'genes5p', 'genes3p',
                          'type', 'distance',
                          'multimap_weighted_encomp_frags',
                          'total_encomp_frags',
                          'unique_alignment_positions',
                          'unique_spanning_frags'])
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
