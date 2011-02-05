'''
Created on Jan 29, 2011

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
import collections
import logging
import os

# local imports
from chimerascan.bx.intersection import Interval, IntervalTree
from chimerascan.lib import config
from chimerascan.lib.feature import GeneFeature
from find_discordant_reads import DiscordantFragment

# constants
MULTIMAP_BINS = (1,2,4,8,16,32,64,128)
CHIMERA_SEP = "|"

class ChimeraMate(object):
    def __init__(self):
        self.start = 0
        self.end = 0
        self.tx_name = None
        self.gene_name = None
        self.strand = "."
        self.exon_start_num = 0
        self.exon_end_num = 0
        self.isize = 0
        self.frac = 0.0
                
class Chimera(object):
    QNAME_COL = 22
    LAST_COL = 24
    SEQ_FIELD_DELIM = ';'

    def __init__(self):
        self.name = None
        self.chimera_type = 0
        self.distance = None 
        self.encompassing_reads = 0
        self.weighted_cov = 0.0
        self.strand_reads = [0, 0]
        self.multimap_cov_hist = None
        self.qnames = []
        self.seqs = []       
        self.mate5p = ChimeraMate()
        self.mate3p = ChimeraMate()
        self.extra_fields = []

    def to_list(self):
        qnames = self.SEQ_FIELD_DELIM.join(x for x in self.qnames)
        seqs1 = self.SEQ_FIELD_DELIM.join(x[0] for x in self.seqs)
        seqs2 = self.SEQ_FIELD_DELIM.join(x[1] for x in self.seqs)
        s = [self.mate5p.tx_name, self.mate5p.start, self.mate5p.end,
             self.mate3p.tx_name, self.mate3p.start, self.mate3p.end,
             self.name, self.weighted_cov, 
             self.mate5p.strand, self.mate3p.strand,
             self.chimera_type, self.distance, self.encompassing_reads, 
             self.strand_reads[0], self.strand_reads[1],
             ','.join(map(str,self.multimap_cov_hist)),
             self.mate5p.isize, self.mate3p.isize,
             '%d-%d' % (self.mate5p.exon_start_num, self.mate5p.exon_end_num),
             '%d-%d' % (self.mate3p.exon_start_num, self.mate3p.exon_end_num),
             self.mate5p.frac, self.mate3p.frac,
             qnames, seqs1, seqs2]
        s.extend(self.extra_fields)
        return s
    
    def from_list(self, fields):
        self.mate5p = ChimeraMate()
        self.mate3p = ChimeraMate()
        self.mate5p.tx_name = fields[0]
        self.mate5p.start = int(fields[1])
        self.mate5p.end = int(fields[2])
        self.mate3p.tx_name = fields[3]
        self.mate3p.start = int(fields[4])
        self.mate3p.end = int(fields[5])
        self.name = fields[6]
        self.mate5p.gene_name = self.name.split(CHIMERA_SEP)[0]
        self.mate3p.gene_name = self.name.split(CHIMERA_SEP)[1]
        self.weighted_cov = float(fields[7])
        self.mate5p.strand = fields[8]
        self.mate3p.strand = fields[9]
        self.chimera_type = fields[10]
        if fields[11] == "None":
            self.distance = None
        else:
            self.distance = int(fields[11])
        self.encompassing_reads = int(fields[12])
        self.strand_reads[0] = int(fields[13])
        self.strand_reads[1] = int(fields[14])        
        self.multimap_cov_hist = map(int,fields[15].split(','))
        self.mate5p.isize = int(fields[16])
        self.mate3p.isize = int(fields[17])        
        self.mate5p.exon_start_num, self.mate5p.exon_end_num = map(int, fields[18].split('-'))
        self.mate3p.exon_start_num, self.mate3p.exon_end_num = map(int, fields[19].split('-'))
        self.mate5p.frac, self.mate3p.frac = map(float, fields[20:22])
        self.qnames = fields[22].split(Chimera.SEQ_FIELD_DELIM)
        self.seqs = zip(fields[23].split(Chimera.SEQ_FIELD_DELIM), 
                        fields[24].split(Chimera.SEQ_FIELD_DELIM))
        self.extra_fields = fields[25:]

    @staticmethod
    def parse(line_iter):
        for line in line_iter:
            fields = line.strip().split('\t')
            c = Chimera()
            c.from_list(fields)
            yield c

def parse_discordant_reads(infh):
    for line in infh:
        fields = line.strip().split('\t')
        frag = DiscordantFragment.from_list(fields)
        yield frag
        
def parse_discordant_gene_pairs(infh):
    prev_tx5p, prev_tx3p = None,None
    reads = []
    for frag in parse_discordant_reads(infh):
#        if frag.discordant_type.is_genome:
#            continue
#        if frag.clust5p.rname == "*":
#            continue
#        if frag.clust3p.rname == "*":
#            continue
        tx5p = frag.clust5p.rname
        tx3p = frag.clust3p.rname
        if (tx5p, tx3p) != (prev_tx5p, prev_tx3p):
            if len(reads) > 0:
                yield prev_tx5p, prev_tx3p, reads
                reads = []
            prev_tx5p, prev_tx3p = tx5p, tx3p
        reads.append(frag)
    if len(reads) > 0:
        yield tx5p, tx3p, reads

def hist(vals, bins):
    '''
    creates a histogram of values in the list of bins provided
    values less than the lowest bin are used in the lowest bin,
    and all values higher than the highest bin are included in the
    highest bin
    '''
    d = collections.defaultdict(lambda: 0)
    for v in vals:
        d[v] += 1
    sorted_keys = sorted(d)
    cov = [0] * len(bins)
    bin_ind = 0
    bin_val = bins[bin_ind]
    count = 0
    for key in sorted_keys:
        while key > bin_val and (bin_ind < len(bins)):
            cov[bin_ind] = count
            count = 0
            bin_ind += 1
            if bin_ind == len(bins):
                break
            bin_val = bins[bin_ind]
        count += d[key]    
    if bin_ind < len(bins):
        cov[bin_ind] += count
    else:
        cov[-1] += count
    return cov

def scoreatpercentile(a, p):
    from math import floor, ceil
    floatind = (len(a)-1) * p
    lowind = int(floor(floatind))
    highind = int(ceil(floatind))
    # interpolate
    val = (1 - (floatind - lowind))*a[lowind] + (1 - (highind - floatind))*a[highind]
    return val

def build_gene_maps(genefile):
    gene_genome_map = {}
    gene_trees = collections.defaultdict(lambda: IntervalTree())    
    # build gene and genome data structures for fast lookup
    for g in GeneFeature.parse(open(genefile)):
        gene_genome_map[g.tx_name] = g
        # add gene to interval tree
        gene_interval = Interval(g.tx_start, g.tx_end, strand=g.strand, value=(g.tx_name))
        gene_trees[g.chrom].insert_interval(gene_interval)
    return gene_genome_map, gene_trees

def get_exon_interval(g, pos):
    exon_iter = reversed(g.exons) if g.strand == '-' else iter(g.exons)
    exon_pos = 0
    exon_num = 0
    for exon_start, exon_end in exon_iter:
        exon_size = exon_end - exon_start
        if exon_pos + exon_size >= pos:
            break
        exon_pos += exon_size
        exon_num += 1    
    if exon_pos + exon_size < pos:
        logging.warning("exon_pos %d + exon_size %d < pos %d - clipping to end of gene" % (exon_pos, exon_size, pos))
    return exon_num, exon_pos, exon_pos + exon_size

def get_chimera_type(fiveprime_gene, threeprime_gene, gene_trees):
    #ORIENTATION_SAME = "Same"
    #ORIENTATION_OPPOSITE = "Opposite"    
    CHIMERA_INTERCHROMOSOMAL = "Interchromosomal"
    CHIMERA_OVERLAP_CONVERGE = "Overlapping_Converging"
    CHIMERA_OVERLAP_DIVERGE = "Overlapping_Diverging"
    CHIMERA_OVERLAP_SAME = "Overlapping_Same"
    CHIMERA_OVERLAP_COMPLEX = "Overlapping_Complex"
    CHIMERA_READTHROUGH = "Read_Through"
    CHIMERA_INTRACHROMOSOMAL = "Intrachromosomal"
    CHIMERA_ADJ_CONVERGE = "Adjacent_Converging"
    CHIMERA_ADJ_DIVERGE = "Adjacent_Diverging"
    CHIMERA_ADJ_SAME = "Adjacent_Same"
    CHIMERA_ADJ_COMPLEX = "Adjacent_Complex"
    #CHIMERA_INTRA_CONVERGE = "Intrachromosomal_Converging"
    #CHIMERA_INTRA_DIVERGE = "Intrachromsomal_Diverging"
    #CHIMERA_INTRA_SAME = "Intrachromosomal_Same"
    CHIMERA_INTRA_COMPLEX = "Intrachromosomal_Complex"
    CHIMERA_UNKNOWN = "Undetermined"
    # get gene information
    chrom1, start5p, end5p, strand1 = fiveprime_gene.chrom, fiveprime_gene.tx_start, fiveprime_gene.tx_end, fiveprime_gene.strand
    chrom2, start3p, end3p, strand2 = threeprime_gene.chrom, threeprime_gene.tx_start, threeprime_gene.tx_end, threeprime_gene.strand
    # interchromosomal
    if chrom1 != chrom2:
        return CHIMERA_INTERCHROMOSOMAL, None    
    # orientation
    same_strand = strand1 == strand2
    # genes on same chromosome so check overlap
    is_overlapping = (start5p < end3p) and (start3p < end5p)            
    if is_overlapping:
        if not same_strand:
            if ((start5p <= start3p and strand1 == "+") or
                (start5p > start3p and strand1 == "-")):                    
                return (CHIMERA_OVERLAP_CONVERGE, 0)
            else:
                return (CHIMERA_OVERLAP_DIVERGE, 0)
        else:
            if ((start5p <= start3p and strand1 == "+") or
                (end5p >= end3p and strand1 == "-")):
                return (CHIMERA_OVERLAP_SAME, 0)
            else:
                return (CHIMERA_OVERLAP_COMPLEX, 0)
    # if code gets here then the genes are on the same chromosome but do not
    # overlap.  first calculate distance (minimum distance between genes)
    if start5p <= start3p:
        distance = start3p - end5p
        between_interval = Interval(end5p, start3p)
    else:
        distance = end3p - start5p
        between_interval = Interval(end3p, start5p)
    # check whether there are genes intervening between the
    # chimera candidates
    genes_between = []
    genes_between_same_strand = []
    for hit in gene_trees[chrom1].find(between_interval.start,
                                       between_interval.end):
        if (hit.start > between_interval.start and
            hit.end < between_interval.end):             
            if hit.strand == strand1:
                genes_between_same_strand.append(hit)
            genes_between.append(hit)

    if same_strand:
        if len(genes_between_same_strand) == 0:
            return CHIMERA_READTHROUGH, distance
        else:
            return CHIMERA_INTRACHROMOSOMAL, distance
    else:
        # check for reads between neighboring genes    
        if len(genes_between) == 0:
            if ((start5p <= start3p and strand1 == "+") or
                (start5p > start3p and strand1 == "-")):                    
                return (CHIMERA_ADJ_CONVERGE, distance)
            elif ((start5p >= start3p and strand1 == "+") or
                  (start5p < start3p and strand1 == "-")):
                return (CHIMERA_ADJ_DIVERGE, distance)
            elif ((start5p <= start3p and strand1 == "+") or
                  (start5p > start3p and strand1 == "-")):
                return (CHIMERA_ADJ_SAME, distance)
            elif ((start5p >= start3p and strand1 == "+") or
                  (start5p < start3p and strand1 == '-')):
                return (CHIMERA_ADJ_COMPLEX, distance)
        else:
            return CHIMERA_INTRA_COMPLEX, distance    
    return CHIMERA_UNKNOWN, distance

#def interval_overlap(chrom1, start1, end1, chrom2, start2, end2):
#    if chrom1 != chrom2:
#        return False
#    return (start1 < end2) and (start2 < end1)    

def uniqueness(r):
    '''return a number representing the multimapping uniqueness of the read'''
    return max(1, r.clust5p.multimaps) * max(1, r.clust3p.multimaps)

def get_chimera_mate(gene, reads, gene_genome_map, trim, is_5prime=True):
    # build 5' and 3' partners
    mate = ChimeraMate()        
    mate.tx_name = gene.tx_name
    # fix gene names with spaces    
    mate.gene_name = '_'.join(gene.gene_name.split())
    starts = []
    ends = []
    for r in reads:
        if is_5prime:
            starts.append(r.clust5p.start)
            ends.append(r.clust5p.end)
        else:
            starts.append(r.clust3p.start)
            ends.append(r.clust3p.end)
    starts = sorted(starts)
    ends = sorted(ends)
    mate.isize = int(scoreatpercentile(ends, 0.95) - scoreatpercentile(starts, 0.05))
    mate.strand = gene.strand
    # TODO: trim to remove unmapped "splash" around exon boundaries
    firststart, lastend = starts[0], ends[-1]
    trimstart = min(firststart + trim, lastend)
    trimend = max(lastend - trim, trimstart + 1)
    # get exon information
    firstexon_num, firstexon_start, firstexon_end = get_exon_interval(gene, trimstart)
    lastexon_num, lastexon_start, lastexon_end = get_exon_interval(gene, trimend)
    if is_5prime:
        mate.start = 0
        mate.end = lastexon_end
    else:
        tx_length = sum(end - start for start, end in gene.exons)
        mate.start = firstexon_start
        mate.end = tx_length
    mate.exon_start_num = firstexon_num
    mate.exon_end_num = lastexon_num
    return mate

def nominate_chimeras(infh, outfh, gene_file, trim=10):
    logging.info("Reading gene information")
    gene_genome_map, gene_trees = build_gene_maps(gene_file)
    logging.info("Processing discordant reads")
    for txname5p, txname3p, reads in parse_discordant_gene_pairs(infh):
        # lookup gene information        
        gene5p = gene_genome_map[txname5p]        
        gene3p = gene_genome_map[txname3p]
        # build 5' and 3' partners
        mate5p = get_chimera_mate(gene5p, reads, gene_genome_map, trim, is_5prime=True)
        mate3p = get_chimera_mate(gene3p, reads, gene_genome_map, trim, is_5prime=False)
        # categorize chimera
        chimera = Chimera()
        chimera.mate5p = mate5p
        chimera.mate3p = mate3p
        chimera_type, distance = get_chimera_type(gene5p, gene3p, gene_trees)
        chimera.chimera_type = chimera_type
        chimera.distance = distance
        # calculate weighted coverage and strand balance
        uniqueness_vals = []
        for r in reads:
            # TODO: assign equal weight to multimapping reads for now,
            # but can use multiple iterations and improve this
            if r.clust1.strand == "+":
                chimera.strand_reads[0] += 1
            else:
                chimera.strand_reads[1] += 1
            uniqueness_vals.append(uniqueness(r))
        chimera.weighted_cov = sum((1.0/x) for x in uniqueness_vals)
        chimera.multimap_cov_hist = hist(uniqueness_vals, MULTIMAP_BINS)
        chimera.encompassing_reads = len(reads)
        # TODO: determine whether fusion likely happens in the intron or 
        # from within the exon itself
        # get qnames and sequences
        for r in reads:
            chimera.qnames.append(r.qname)
            chimera.seqs.append((r.clust1.seq, r.clust2.seq))
        chimera.name = '%s%s%s' % (mate5p.gene_name, CHIMERA_SEP, mate3p.gene_name)
        fields = chimera.to_list()
        print >>outfh, '\t'.join(map(str, fields))

def main():
    from optparse import OptionParser
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = OptionParser("usage: %prog [options] <sortedchimeras.bedpe> <chimeras.txt>")
    parser.add_option("--index", dest="index_dir",
                      help="Path to chimerascan index directory")
    options, args = parser.parse_args()
    input_file = args[0]
    output_file = args[1]
    gene_file = os.path.join(options.index_dir, config.GENE_FEATURE_FILE)
    nominate_chimeras(open(input_file), 
                      open(output_file,'w'),
                      gene_file)


if __name__ == '__main__':
    main()
