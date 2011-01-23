'''
Created on Jan 22, 2011

@author: mkiyer
'''
'''
Created on Jan 11, 2011

@author: mkiyer
'''
import collections
import itertools
import logging
import operator
import os
import sys

# local libs
import pysam
from bx.intersection import Interval, IntervalTree
from bx.cluster import ClusterTree

# local imports
import config
from feature import GeneFeature
from seq import DNA_reverse_complement
from gene_to_genome import build_gene_maps, get_gene_tids
from alignment_parser import parse_pe_sam_file

def build_exon_trees(samfh, genefile):
    rname_tid_map = dict((rname,i) for i,rname in enumerate(samfh.references))
    exon_trees = collections.defaultdict(lambda: IntervalTree())    
    # build gene and genome data structures for fast lookup
    for g in GeneFeature.parse(open(genefile)):
        name = config.GENE_REF_PREFIX + g.tx_name
        if name not in rname_tid_map:
            continue
        if g.chrom not in rname_tid_map:
            continue
        gene_tid = rname_tid_map[name]
        # get reference index in sam file
        chrom_tid = rname_tid_map[g.chrom]        
        # add gene to interval tree
        for start,end in g.exons[1::-1]:        
            exon_interval = Interval(start, end, chrom=chrom_tid, strand=g.strand, value=gene_tid)
            exon_trees[chrom_tid].insert_interval(exon_interval)
    return dict(exon_trees)


def realign_read(read, exon_trees):    
    for hit in exon_trees[read.rname].find(read.pos, read.aend):
        print hit

def realign_split_reads(reads, gene_tid_list, exon_trees):
    for r in reads: 
        if gene_tid_list[r.rname] is None:
            realign_read(r, exon_trees)
        yield r

def realign_genome_reads(input_bam_file, output_bam_file, gene_file):
    # build a map of gene name to genome coords
    logging.info("Reading gene index")
    infh = pysam.Samfile(input_bam_file, "rb")    
    gene_tid_list = get_gene_tids(infh)
    exon_trees = build_exon_trees(infh, gene_file)
    outfh = pysam.Samfile("-", "w", template=infh)
    #outfh = pysam.Samfile(output_bam_file, "wb", template=infh)
    for pe_reads in parse_pe_sam_file(infh):
        for mate_partitions in pe_reads:
            for splits in mate_partitions:
                for reads in splits:
                    for r in realign_split_reads(reads, gene_tid_list, exon_trees):
                        outfh.write(r)


def main():
    from optparse import OptionParser
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = OptionParser("usage: %prog [options] <bam> <out.bedpe>")
    parser.add_option("--index", dest="index_dir",
                      help="Path to chimerascan index directory")
    options, args = parser.parse_args()
    if options.index_dir is None:
        parser.error("Must specify chimerascan index directory with --index")
    input_bam_file = args[0]
    output_bam_file = args[1]    
    gene_file = os.path.join(options.index_dir, config.GENE_FEATURE_FILE)    
    realign_genome_reads(input_bam_file, output_bam_file, gene_file)

if __name__ == '__main__':
    main()
