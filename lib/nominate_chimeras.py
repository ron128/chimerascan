'''
Created on Jan 17, 2011

@author: mkiyer
'''
import collections
import logging
import os

# local imports
from bx.intersection import Interval, IntervalTree
import config
from feature import GeneFeature
from find_discordant_reads import Chimera

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
    assert exon_pos + exon_size >= pos
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

def interval_overlap(chrom1, start1, end1, chrom2, start2, end2):
    if chrom1 != chrom2:
        return False
    return (start1 < end2) and (start2 < end1)    

def parse_discordant_reads(infh):
    prev_tx5p, prev_tx3p = None,None
    reads = []
    for line in infh:
        chimera = Chimera.from_bedpe(line)
        tx5p = chimera.mate5p.interval.chrom
        tx3p = chimera.mate3p.interval.chrom       
        if (tx5p, tx3p) != (prev_tx5p, prev_tx3p):
            if len(reads) > 0:
                yield prev_tx5p, prev_tx3p, reads
                reads = []
            prev_tx5p, prev_tx3p = tx5p, tx3p
        reads.append(chimera)
    if len(reads) > 0:
        yield tx5p, tx3p, reads

def nominate_chimeras(infh, outfh, gene_file,
                      contam_refs, trim=10):
    logging.info("Reading gene information")
    gene_genome_map, gene_trees = build_gene_maps(gene_file)
    logging.info("Processing discordant reads")
    for txname5p, txname3p, reads in parse_discordant_reads(infh):
        # get start/end coordinates of chimeric reads
        start5p = min(r.mate5p.interval.start for r in reads)
        end5p = max(r.mate5p.interval.end for r in reads) 
        start3p = min(r.mate3p.interval.start for r in reads)
        end3p = max(r.mate3p.interval.end for r in reads)
        # TODO: trim to allow unmapped "splash" around exon boundaries
        trimstart5p = min(start5p + trim, end5p)
        trimend5p = max(end5p - trim, start5p) 
        trimstart3p = min(start3p + trim, end3p)
        trimend3p = max(end3p - trim, start3p)        
        # lookup gene information
        gene5p = gene_genome_map[txname5p]
        gene3p = gene_genome_map[txname3p]
        gene_length5p = sum(end - start for start, end in gene5p.exons)
        gene_length3p = sum(end - start for start, end in gene3p.exons)                
        # categorize chimera
        chimera_type, distance = get_chimera_type(gene5p, gene3p, gene_trees)
        # get exon information
        firstexon5p_num, firstexon5p_start, firstexon5p_end = get_exon_interval(gene5p, trimstart5p)
        lastexon5p_num, lastexon5p_start, lastexon5p_end = get_exon_interval(gene5p, trimend5p)
        firstexon3p_num, firstexon3p_start, firstexon3p_end = get_exon_interval(gene3p, trimstart3p)
        lastexon3p_num, lastexon3p_start, lastexon3p_end = get_exon_interval(gene3p, trimend3p)
        # TODO: determine whether fusion likely happens in the intron or 
        # from within the exon itself
        # get qnames and sequences
        qnames = ';'.join(r.qname for r in reads)
        seqs5p = ';'.join(r.mate5p.seq for r in reads)
        seqs3p = ';'.join(r.mate3p.seq for r in reads)
        # fix gene names with spaces
        name5p = '_'.join(gene5p.gene_name.split())
        name3p = '_'.join(gene3p.gene_name.split())
        s = [gene5p.tx_name, lastexon5p_end, gene_length5p,
             gene3p.tx_name, firstexon3p_start, gene_length3p,
             '%s-%s' % (name5p, name3p),
             len(reads),
             gene5p.strand, gene3p.strand,
             chimera_type, distance,
             '%d-%d' % (firstexon5p_num, lastexon5p_num),
             '%d-%d' % (firstexon3p_num, lastexon3p_num),
             qnames, seqs5p, seqs3p]
        print >>outfh, '\t'.join(map(str, s))

def main():
    from optparse import OptionParser
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = OptionParser("usage: %prog [options] <sortedchimeras.bedpe> <chimeras.txt>")
    parser.add_option("--index", dest="index_dir",
                      help="Path to chimerascan index directory")
    parser.add_option("--contam-refs", dest="contam_refs", default=None)
    options, args = parser.parse_args()
    input_file = args[0]
    output_file = args[1]
    gene_file = os.path.join(options.index_dir, config.GENE_FEATURE_FILE)
    nominate_chimeras(open(input_file), open(output_file,'w'), 
                      gene_file, options.contam_refs)

if __name__ == '__main__':
    main()
