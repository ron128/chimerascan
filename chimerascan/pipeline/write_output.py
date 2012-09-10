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
import argparse
import logging
import os
import sys
import collections
import shelve

from chimerascan.bx.intersection import Interval, IntervalTree

from chimerascan.lib import config
from chimerascan.lib.chimera import Chimera, \
    parse_discordant_cluster_pair_file, get_chimera_type
from chimerascan.lib.feature import TranscriptFeature

def build_genome_transcript_trees(transcripts):
    genome_tx_trees = collections.defaultdict(lambda: IntervalTree())    
    transcript_dict = {}
    for t in transcripts:
        # add to dict
        transcript_dict[t.tx_id] = t
        # add exons to interval tree
        for start,end in t.exons:
            interval = Interval(start, end, strand=t.strand, value=t.tx_id)
            genome_tx_trees[t.chrom].insert_interval(interval)
    return transcript_dict, genome_tx_trees

def lookup_transcripts(cluster, transcript_dict, genome_tx_trees):
    tx_ids = set()
    for hit in genome_tx_trees[cluster.rname].find(cluster.start, cluster.end):
        if hit.strand != cluster.strand:
            continue        
        tx_ids.add(hit.value)
    hits = [transcript_dict[tx_id] for tx_id in sorted(tx_ids)]
    return hits

def get_transcript_info(transcripts, annotation_source):
    tx_names = set() 
    gene_names = set()
    biotypes = set()
    for t in transcripts:
        if annotation_source in t.annotation_sources:
            i = t.annotation_sources.index(annotation_source)
            tx_names.add(t.tx_names[i])
            gene_names.add(t.gene_names[i])
        if t.gene_biotype != 'na':
            biotypes.add(t.gene_biotype)
    if len(tx_names) == 0:
        tx_names.add('na')
    if len(gene_names) == 0:
        gene_names.add('na')
    if len(biotypes) == 0:
        biotypes.add('na')
    return tx_names, gene_names, biotypes

def make_chimera(cluster_pair, 
                 cluster_shelve,
                 transcript_dict,
                 genome_tx_trees,
                 annotation_source):
    # lookup 5' and 3' clusters
    cluster5p = cluster_shelve[str(cluster_pair.id5p)]
    cluster3p = cluster_shelve[str(cluster_pair.id3p)]
    # get 5' and 3' transcripts
    transcripts5p = lookup_transcripts(cluster5p, transcript_dict, genome_tx_trees)
    transcripts3p = lookup_transcripts(cluster3p, transcript_dict, genome_tx_trees)
    # lookup chimera type and distance
    chimera_type, distance = get_chimera_type(cluster5p, cluster3p, 
                                              transcripts5p, transcripts3p, 
                                              transcript_dict, genome_tx_trees)
    # format transcript information
    tx_names_5p, gene_names_5p, biotypes_5p = get_transcript_info(transcripts5p, annotation_source)
    tx_names_3p, gene_names_3p, biotypes_3p = get_transcript_info(transcripts3p, annotation_source)
    # make chimera object
    c = Chimera()
    c.rname5p = cluster5p.rname
    c.start5p = cluster5p.start
    c.end5p = cluster5p.end
    c.rname3p = cluster3p.rname
    c.start3p = cluster3p.start
    c.end3p = cluster3p.end
    c.chimera_id = "CHIMERA%d" % (cluster_pair.pair_id)
    frags = set(cluster_pair.qnames)
    frags.update(cluster_pair.spanning_qnames)
    c.num_frags = len(frags)
    c.strand5p = cluster5p.strand
    c.strand3p = cluster3p.strand
    c.chimera_type = chimera_type
    c.distance = distance
    c.num_discordant_frags = len(cluster_pair.qnames)
    c.num_spanning_frags = len(cluster_pair.spanning_qnames)
    c.num_discordant_frags_5p = len(cluster5p.qnames)
    c.num_discordant_frags_3p = len(cluster3p.qnames)
    c.num_concordant_frags_5p = cluster5p.concordant_frags
    c.num_concordant_frags_3p = cluster3p.concordant_frags
    c.biotypes_5p = sorted(biotypes_5p)
    c.biotypes_3p = sorted(biotypes_3p)
    c.genes_5p = sorted(gene_names_5p)
    c.genes_3p = sorted(gene_names_3p)
    c.transcripts_5p = sorted(tx_names_5p)
    c.transcripts_3p = sorted(tx_names_3p)
    return c

def write_output(transcripts, cluster_shelve_file, cluster_pair_file, 
                 read_name_file, output_file, 
                 annotation_source="ensembl"):
    # load cluster and read name database files
    cluster_shelve = shelve.open(cluster_shelve_file, 'r')
    read_name_fh = open(read_name_file, 'r')   
    # map genome coordinates to transcripts
    logging.debug("Creating mapping between genome coordinates and transcripts")
    transcript_dict, genome_tx_trees = build_genome_transcript_trees(transcripts)
    logging.debug("Writing output")
    outfh = open(output_file, "w")
    print >>outfh, '#' + '\t'.join(Chimera._fields)
    for cluster_pair in parse_discordant_cluster_pair_file(open(cluster_pair_file)):
        c = make_chimera(cluster_pair, cluster_shelve, transcript_dict, 
                         genome_tx_trees, annotation_source)
        print >>outfh, str(c)
    # cleanup
    outfh.close()
    read_name_fh.close()
    cluster_shelve.close()
    return config.JOB_SUCCESS

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("--ann", dest="annotation_source", default="ensembl")
    parser.add_argument("transcript_file")
    parser.add_argument("cluster_shelve_file")
    parser.add_argument("cluster_pair_file")
    parser.add_argument("read_name_file")
    parser.add_argument("output_file")
    args = parser.parse_args()    
    # read transcript features
    logging.debug("Reading transcript features")
    transcripts = list(TranscriptFeature.parse(open(args.transcript_file)))
    # run main function
    retcode = write_output(transcripts, args.cluster_shelve_file, 
                           args.cluster_pair_file, args.read_name_file, 
                           args.output_file, args.annotation_source)
    return retcode

if __name__ == "__main__":
    sys.exit(main())
