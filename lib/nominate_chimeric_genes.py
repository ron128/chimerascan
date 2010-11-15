'''
Created on Nov 7, 2010

@author: mkiyer
'''
import logging
import argparse
import collections
import operator
import itertools

import pysam
from bx.intervals.cluster import ClusterTree
from bx.intervals.intersection import Interval, IntervalTree

from base import parse_multihit_sam_file, get_aligned_read_intervals

class BEDFeature(object):
    __slots__ = ('chrom', 'tx_start', 'tx_end', 'name', 'score', 
                 'strand', 'genomic_exons', 'tx_exons', 'tx_length',
                 'alias')
    
    def __str__(self):
        return '\t'.join(map(str, [self.chrom, self.tx_start, self.tx_end, 
                                   self.name, self.score, self.strand,
                                   self.genomic_exons, self.tx_exons,
                                   self.tx_length]))
    
    @staticmethod
    def from_bed12_string(line):
        if line is None:
            return None
        line = line.strip()
        if line.startswith('#'):
            logging.debug("skipping comment line: %s" % (line))
            return None
        if line.startswith('track'):
            logging.debug("skipping track header line: %s"  % (line))
            return None
        fields = line.split('\t')
        # first six fields are required
        g = BEDFeature()
        g.chrom = fields[0]
        g.tx_start = int(fields[1])
        g.tx_end = int(fields[2])
        g.name = fields[3]
        g.score = fields[4]
        g.strand = fields[5]
        exon_sizes = map(int, fields[10].split(',')[:-1])
        exon_starts = map(int, fields[11].split(',')[:-1])        
        g.genomic_exons = []
        for e_start, e_size in itertools.izip(exon_starts, exon_sizes):
            genomic_e_start = g.tx_start + e_start
            g.genomic_exons.append((genomic_e_start, genomic_e_start + e_size))
        # account for strand
        if g.strand == "-":
            exon_sizes.reverse()
            g.genomic_exons.reverse()
        tx_pos = 0
        g.tx_length = sum(exon_sizes)
        g.tx_exons = []
        for e_size in exon_sizes:
            g.tx_exons.append((tx_pos, tx_pos + e_size))
            tx_pos += e_size
        if g.strand == "-":
            assert g.genomic_exons[0][1] == g.tx_end
            assert g.genomic_exons[-1][0] == g.tx_start
        else:
            assert g.genomic_exons[0][0] == g.tx_start
            assert g.genomic_exons[-1][1] == g.tx_end            
        assert g.tx_exons[-1][1] == g.tx_length
        return g

def read_gene_alias_file(alias_file):
    alias_lookup = {}
    for line in open(alias_file):
        if not line or not line.strip():
            continue
        if line.startswith('#'):
            continue
        fields = line.strip().split('\t')
        alias_lookup[fields[0]] = fields[1]
    return alias_lookup

def build_genome_gene_map(samfh, bedfile, alias_file=None):
    rname_tid_map = dict((rname,i) for i,rname in enumerate(samfh.references))    
    trees = collections.defaultdict(lambda: IntervalTree())    
    intervals = collections.defaultdict(lambda: [])
    exon_gene_map = {}
    current_exon_id = 0

    if alias_file is not None:
        alias_lookup = read_gene_alias_file(alias_file)
    else:
        alias_lookup = {}

    for line in open(bedfile):
        g = BEDFeature.from_bed12_string(line)
        if g.chrom not in rname_tid_map:
            continue
        g.alias = alias_lookup.get(g.name, g.name)
        tid = rname_tid_map[g.chrom]
        for exon_num in xrange(len(g.genomic_exons)):
            genomic_e_start, genomic_e_end = g.genomic_exons[exon_num]
            exon_key = (g.chrom, genomic_e_start, genomic_e_end, g.strand)            
            if exon_key in exon_gene_map:
                exon_id = exon_gene_map[exon_key]
            else:
                exon_id = current_exon_id
                current_exon_id += 1
                exon_gene_map[exon_key] = exon_id
                trees[tid].insert_interval(Interval(genomic_e_start, genomic_e_end, 
                                                    strand=g.strand, chrom=g.chrom, value=exon_id))
            intervals[exon_id].append((g, exon_num))
            #print 'exon id', exon_id, 'genes', [(x[0].name,x[0].genomic_exons[x[1]]) for x in intervals[exon_id]]
    return trees, intervals

def intersect_read_multihits(reads, exon_trees, contam_tids):
    # find genes overlapping the reads
    hits = [collections.defaultdict(lambda: []),
            collections.defaultdict(lambda: [])]
    for read in reads:
        if read.is_unmapped:
            continue
        if read.rname in contam_tids:
            continue
        read_strand = "-" if read.is_reverse else "+"
        mate = 0 if read.is_read1 else 1
        for interval in get_aligned_read_intervals(read):
            for hit in exon_trees[read.rname].find(interval[0], interval[1]):
                # ensure exon entirely encompasses the read
                if hit.start <= interval[0] and hit.end >= interval[1]:
                    hits[mate][(hit.value, hit.strand)].append(read_strand)
    return hits

def cluster_reads(samfh, exon_trees, exon_intervals, contam_tids, min_dist):    
    debug_every = 1e5
    debug_next = debug_every    
    read_num = 0    
    exon_pair_reads = collections.defaultdict(lambda: set())
    try:
        sam_iter = parse_multihit_sam_file(samfh)
        while True:
            reads = sam_iter.next()            
            # logging debug output
            read_num += 1
            if read_num == debug_next:
                debug_next += debug_every
                logging.debug("Finished clustering %d reads" % read_num)            
            # find genes overlapping the reads
            mate1_hits, mate2_hits = intersect_read_multihits(reads, exon_trees, contam_tids)
            # count cluster-pair hits
            for exon_data1,read_strands1 in mate1_hits.iteritems():
                exon_id1, exon_strand1 = exon_data1                
                for exon_data2,read_strands2 in mate2_hits.iteritems():
                    exon_id2, exon_strand2 = exon_data2                                    
                    # skip when map to the same gene
                    if exon_id1 == exon_id2:
                        continue                    
                    # ensure read strandedness matches exon strandedness
                    valid_strands = False
                    for read_strand1 in read_strands1:
                        for read_strand2 in read_strands2:
                            if (((exon_strand1 == exon_strand2) and (read_strand1 != read_strand2)) or
                                ((exon_strand1 != exon_strand2) and (read_strand1 == read_strand2))):
                                valid_strands = True
                                break
                        if valid_strands:
                            break
                    if not valid_strands:
                        continue
                    # figure out 5'/3' orientation based on concordance between
                    # read strand and exon strand
                    if read_strand1 == exon_strand1:
                        exon_pair = (exon_id1, exon_id2)
                    else:
                        exon_pair = (exon_id2, exon_id1)
                    exon_pair_reads[exon_pair].add(reads[0].qname)
    except StopIteration:
        pass
    return exon_pair_reads

def filter_chimeric_genes(exon_pair_reads, exon_intervals):
    for exon_ids,read_names in exon_pair_reads.iteritems():
        a_genes = exon_intervals[exon_ids[0]]
        b_genes = exon_intervals[exon_ids[1]]
        for a, a_exon_num in a_genes:
            for b, b_exon_num in b_genes:
                print '\t'.join(map(str, [a.name, a.tx_exons[a_exon_num][1], a.tx_length,
                                          b.name, b.tx_exons[b_exon_num][0], b.tx_length,
                                          '-'.join([a.alias, b.alias]), len(read_names),
                                          a.strand, b.strand, ','.join(read_names)]))
#                print '\t'.join(map(str, [b.name, b.tx_exons[b_exon_num][1], b.tx_length,
#                                          a.name, a.tx_exons[a_exon_num][0], a.tx_length,
#                                          '-'.join([b.alias, a.alias]), num_reads,
#                                          b.strand, a.strand, "read_names"]))

#                print '\t'.join(map(str, [b.name, b.tx_exons[b_exon_num][1], b.tx_length,
#                                          a.name, a.tx_exons[a_exon_num][0], a.tx_length,
#                                          "hugo", b.strand, a.strand, num_reads, "read_names"]))
#                print '\t'.join(map(str, [a.chrom, a.start, a.end, cluster_ids[0], num_reads, a.strand, 
#                                          b.chrom, b.start, b.end, cluster_ids[1], num_reads, b.strand]))
#                if a.chrom == b.chrom:
#                    dist = min(abs(a.start - b.start), abs(a.end - b.start),
#                               abs(a.start - b.end), abs(a.end - b.end))
#                    if dist < min_dist:
#                        continue

def get_tids(samfh, rnames):
    rname_tid_map = dict((rname,i) for i,rname in enumerate(samfh.references))    
    return [rname_tid_map[rname] for rname in rnames]

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")    
    parser = argparse.ArgumentParser()
    parser.add_argument("--min-dist", dest="min_dist", type=int, default=300000)
    parser.add_argument("--library-type", dest="library_type", default="fr")
    parser.add_argument("--aliases", dest="alias_file", default=None)
    parser.add_argument("--contam", dest="contam_refs", default=None)
    parser.add_argument("bed_file")
    parser.add_argument("sam_file")
    options = parser.parse_args()

    bamfh = pysam.Samfile(options.sam_file, "rb")
    if options.contam_refs is None:
        contam_tids = set()
    else:
        contam_tids = set(get_tids(bamfh, options.contam_refs.split(',')))
    # build interval trees from annotated exons
    logging.debug("Building exon intersection tables")
    exon_trees, exon_intervals = build_genome_gene_map(bamfh, options.bed_file, options.alias_file)
    logging.debug("Clustering chimeric gene pairs")
    exon_pair_reads = cluster_reads(bamfh, exon_trees, exon_intervals, contam_tids, options.min_dist)    
    bamfh.close()
    logging.debug("Filtering chimeric genes")
    filter_chimeric_genes(exon_pair_reads, exon_intervals)    

if __name__ == '__main__': main()

#            for r in reads:
#                print r.qname, samfh.getrname(r.rname) if r.rname >= 0 else "NA", r.rname, r.pos, r.is_read1, r.is_read2
#            print r1_hits, r2_hits
#            print '-------'
#    sam_iter = parse_multihit_sam_file(bamfh)
#    for reads in sam_iter:
#        for r in reads:
#            print r.qname
#        print '-----'
#    return