'''
Created on Jan 22, 2011

@author: mkiyer
'''
'''
Created on Jan 11, 2011

@author: mkiyer
'''
import collections
import logging
import os

# local libs
from chimerascan import pysam
from chimerascan.bx.intersection import Interval, IntervalTree
from chimerascan.bx.cluster import ClusterTree
from chimerascan.lib import config
from chimerascan.lib.feature import GeneFeature
from chimerascan.lib.seq import DNA_reverse_complement
from chimerascan.lib.gene_to_genome import build_gene_maps, get_gene_tids
from chimerascan.lib.alignment_parser import parse_pe_sam_file

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
#def rescue_genome_mappings(pe_reads, gene_tid_list, gene_trees):
#    for mate, partitions in enumerate(pe_reads):
#        for partition in partitions:
#            splits5p = []
#            splits3p = []
#            codes5p = set()
#            codes3p = set()
#            for split_reads in partition:
#                reads5p = []
#                reads3p = []
#                # TODO: we select reads in the best 'strata', that is, the
#                # set of reads with fewest mismatches to the reference.  is
#                # this the best strategy, or should other metrics be employed
#                # to choose from among multimapping reads?
#                best_reads = select_best_mismatches(split_reads)
#                for r in best_reads:
#                    if r.is_unmapped:
#                        reads3p.append(r)
#                        reads5p.append(r)
#                        if r.rname != -1:
#                            code = GENOME
#                        elif r.opt('XM') > 0:
#                            code = MULTIMAP
#                        else:
#                            code = NM
#                        codes5p.add(code)
#                        codes3p.add(code)
#                    else:
#                        if r.is_reverse:
#                            # determine sense/antisense by assuming that
#                            # 5' reads are sense and 3' reads are antisense                    
#                            reads3p.append(r)
#                            codes3p.add(MAP)
#                        else:
#                            reads5p.append(r)
#                            codes5p.add(MAP)
#                splits5p.append(reads5p)
#                splits3p.append(reads3p)
#            if all(len(reads) > 0 for reads in splits5p):
#                yield 0, codes5p, splits5p
#            if all(len(reads) > 0 for reads in splits3p):
#                yield 1, codes3p, splits3p    
#            # TODO: for now, we treat genomic reads as unmapped because they
#            # require further processing.  by treating unmapped we allow them
#            # to be consider as mis-mapped spanning reads in future steps,
#            # and allow the other segments in the paired-end fragment to 
#            # determine  
#            if tid_list[read.rname] is None:
#                read.is_unmapped = True

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
