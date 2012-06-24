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

import pysam

from chimerascan.lib.chimera import Chimera, get_chimera_type
from chimerascan.lib.feature import TranscriptFeature
from chimerascan.lib import config
from chimerascan.lib.transcriptome import build_transcript_genome_map, \
    build_genome_transcript_trees, build_transcript_map, transcript_to_genome_pos
from chimerascan.pipeline.filter_chimeras import get_wildtype_frags

def get_chimera_groups(input_file, tx_id_map):
    # group chimeras in the same genomic cluster with the same
    # breakpoint
    cluster_chimera_dict = collections.defaultdict(lambda: [])
    for c in Chimera.parse(open(input_file)):
        # get cluster of overlapping genes
        cluster5p = tx_id_map[c.tx_name_5p].cluster_id
        cluster3p = tx_id_map[c.tx_name_3p].cluster_id
        # add to dictionary
        cluster_chimera_dict[(cluster5p,cluster3p)].append(c)
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

def write_output(input_file, bam_file, output_file, index_dir):
    # read transcripts
    logging.debug("Reading transcripts")
    transcript_file = os.path.join(index_dir, config.TRANSCRIPT_FEATURE_FILE)
    transcripts = list(TranscriptFeature.parse(open(transcript_file)))
    # build a lookup table to get genome coordinates from transcript 
    # coordinates
    transcript_genome_map = build_transcript_genome_map(transcripts)
    tx_id_map = build_transcript_map(transcripts)
    genome_tx_trees = build_genome_transcript_trees(transcripts)
    # open BAM file for checking wild-type isoform
    bamfh = pysam.Samfile(bam_file, "rb")   
    # group chimera isoforms together
    lines = []
    chimera_clusters = 0
    for key,chimeras in get_chimera_groups(input_file, tx_id_map):
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
        chimera_type, distance = get_chimera_type(tx_id_map[c.tx_name_5p],
                                                  tx_id_map[c.tx_name_3p],
                                                  genome_tx_trees)
        # get genomic positions of chimera
        chrom5p,strand5p,start5p = transcript_to_genome_pos(c.tx_name_5p, c.tx_start_5p, transcript_genome_map)
        chrom5p,strand5p,end5p = transcript_to_genome_pos(c.tx_name_5p, c.tx_end_5p-1, transcript_genome_map)
        if strand5p == 1:
            start5p,end5p = end5p,start5p
        chrom3p,strand3p,start3p = transcript_to_genome_pos(c.tx_name_3p, c.tx_start_3p, transcript_genome_map)
        chrom3p,strand3p,end3p = transcript_to_genome_pos(c.tx_name_3p, c.tx_end_3p-1, transcript_genome_map)
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
        # get isoform fraction
        num_wt_frags_5p, num_wt_frags_3p = get_wildtype_frags(c, bamfh)
        num_chimeric_frags = c.get_num_frags()
        frac5p = float(num_chimeric_frags) / (num_chimeric_frags + num_wt_frags_5p)
        frac3p = float(num_chimeric_frags) / (num_chimeric_frags + num_wt_frags_3p)
        # setup fields of BEDPE file
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
                  frac5p, frac3p,
                  ','.join(spanning_fasta_lines),
                  ','.join(names)]
        lines.append(fields)
        chimera_clusters += 1
    bamfh.close()
    logging.debug("Clustered chimeras: %d" % (chimera_clusters))
    # sort
    lines = sorted(lines, key=operator.itemgetter(18, 17, 16), reverse=True)    
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
                          'isoform_fraction_5p',
                          'isoform_fraction_3p',
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
    parser = OptionParser("usage: %prog [options] <index_dir> <in.txt> <bam_file> <out.txt>")
    options, args = parser.parse_args()
    index_dir = args[0]
    input_file = args[1]
    bam_file = args[2]
    output_file = args[3]
    return write_output(input_file, bam_file, output_file, index_dir)

if __name__ == "__main__":
    sys.exit(main())
