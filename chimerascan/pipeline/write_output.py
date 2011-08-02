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

from chimerascan.lib.chimera import Chimera, get_chimera_type
from chimerascan.lib import config
from chimerascan.lib.gene_to_genome import build_gene_to_genome_map, \
    build_tx_cluster_map, build_genome_tx_trees, build_tx_name_gene_map, gene_to_genome_pos

def get_chimera_groups(input_file, gene_file):
    # build a lookup table to get gene clusters from transcript name    
    tx_cluster_map = build_tx_cluster_map(open(gene_file))
    # build a lookup table to get genome coordinates from transcript 
    # coordinates
    # TODO: can either group by exact breakpoint, or just by
    # gene cluster
    #tx_genome_map = build_gene_to_genome_map(open(gene_file))
    # group chimeras in the same genomic cluster with the same
    # breakpoint
    cluster_chimera_dict = collections.defaultdict(lambda: [])
    for c in Chimera.parse(open(input_file)):
        # get cluster of overlapping genes
        cluster5p = tx_cluster_map[c.tx_name_5p]
        cluster3p = tx_cluster_map[c.tx_name_3p]
        # get genomic positions of breakpoints
        #coord5p = gene_to_genome_pos(c.partner5p.tx_name, c.partner5p.end-1, tx_genome_map)
        #coord3p = gene_to_genome_pos(c.partner3p.tx_name, c.partner3p.start, tx_genome_map)
        # add to dictionary
        cluster_chimera_dict[(cluster5p,cluster3p)].append(c)
        # TODO: use this grouping instead?
        #cluster_chimera_dict[(cluster5p,cluster3p,coord5p,coord3p)].append(c)
    for key,chimeras in cluster_chimera_dict.iteritems():
        yield key,chimeras

def get_best_coverage_chimera(chimeras):
    stats = []
    for c in chimeras:
        # TODO: come up with a way to prioritize here (spanning included?)
        stats.append((c,
                      c.get_num_unique_positions(),
                      c.get_num_frags()))
    sorted_stats = sorted(stats, key=operator.itemgetter(1,2), reverse=True)
    return sorted_stats[0][0]

def write_output(input_file, output_file, index_dir):
    gene_file = os.path.join(index_dir, config.GENE_FEATURE_FILE)
    # build a lookup table to get genome coordinates from transcript 
    # coordinates
    tx_genome_map = build_gene_to_genome_map(open(gene_file))    
    tx_name_gene_map = build_tx_name_gene_map(gene_file)    
    genome_tx_trees = build_genome_tx_trees(gene_file)
    # group chimera isoforms together
    # TODO: requires reading all chimeras into memory
    lines = []
    chimera_clusters = 0
    for key,chimeras in get_chimera_groups(input_file, gene_file):
        #cluster5p,cluster3p,coord5p,coord3p = key
        #cluster5p,cluster3p = key
        txs5p = set()
        txs3p = set()
        genes5p = set()
        genes3p = set()
        names = set()
        for c in chimeras:
            txs5p.add("%s:%d-%d" % (c.tx_name_5p, c.tx_start_5p, c.tx_end_5p-1))
            txs3p.add("%s:%d-%d" % (c.tx_name_3p, c.tx_start_3p, c.tx_end_3p-1))
            genes5p.add(c.gene_name_5p)
            genes3p.add(c.gene_name_3p)
            names.add(c.name)
        c = get_best_coverage_chimera(chimeras)
        # get chimera type and distance between genes
        chimera_type, distance = get_chimera_type(tx_name_gene_map[c.tx_name_5p],
                                                  tx_name_gene_map[c.tx_name_3p],
                                                  genome_tx_trees)
        # get genomic positions of chimera
        chrom5p,strand5p,start5p = gene_to_genome_pos(c.tx_name_5p, c.tx_start_5p, tx_genome_map)
        chrom5p,strand5p,end5p = gene_to_genome_pos(c.tx_name_5p, c.tx_end_5p-1, tx_genome_map)
        if strand5p == 1:
            start5p,end5p = end5p,start5p
        chrom3p,strand3p,start3p = gene_to_genome_pos(c.tx_name_3p, c.tx_start_3p, tx_genome_map)
        chrom3p,strand3p,end3p = gene_to_genome_pos(c.tx_name_3p, c.tx_end_3p-1, tx_genome_map)
        if strand3p == 1:
            start3p,end3p = end3p,start3p
        # get breakpoint spanning sequences
        spanning_seqs = set()
        spanning_fasta_lines = []
        for dr in c.get_spanning_reads():
            if dr.seq in spanning_seqs:
                continue
            spanning_seqs.add(dr.seq)
            spanning_fasta_lines.extend([">%s/%d;pos=%d;strand=%s" % 
                                         (dr.qname, dr.readnum+1, dr.pos, 
                                          "-" if dr.is_reverse else "+"), 
                                         dr.seq])           
        fields = [chrom5p, start5p, end5p,
                  chrom3p, start3p, end3p,
                  "CLUSTER%d" % (chimera_clusters),
                  c.get_num_frags(),
                  "+" if (strand5p == 0) else "-",
                  "+" if (strand3p == 0) else "-",
                  ','.join(txs5p),
                  ','.join(txs3p),
                  ','.join(genes5p),
                  ','.join(genes3p),
                  chimera_type, distance,
                  c.get_num_frags(),
                  c.get_num_spanning_frags(),
                  c.get_num_unique_positions(),
                  ','.join(spanning_fasta_lines),
                  ','.join(names)]
        lines.append(fields)
        chimera_clusters += 1
    logging.debug("Clustered chimeras: %d" % (chimera_clusters))
    # sort
    lines = sorted(lines, key=operator.itemgetter(16, 12, 13), reverse=True)    
    f = open(output_file, "w")
    print >>f, '\t'.join(['#chrom5p', 'start5p', 'end5p', 
                          'chrom3p', 'start3p', 'end3p',
                          'chimera_cluster_id', 'score', 
                          'strand5p', 'strand3p',
                          'transcript_ids_5p', 'transcript_ids_3p',
                          'genes5p', 'genes3p',
                          'type', 'distance',
                          'total_frags', 
                          'spanning_frags',
                          'unique_alignment_positions',
                          'breakpoint_spanning_reads',
                          'chimera_ids'])
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
