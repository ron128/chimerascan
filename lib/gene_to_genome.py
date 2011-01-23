'''
Created on Oct 25, 2010

@author: mkiyer
'''
import collections
import logging

import cPickle as pickle

# local libraries
from bx.intersection import Interval, IntervalTree
import pysam

# local imports
import config
from seq import DNA_reverse_complement
from feature import GeneFeature
from base import parse_multihit_alignments

def build_gene_maps(samfh, genefile):
    rname_tid_map = dict((rname,i) for i,rname in enumerate(samfh.references))
    gene_genome_map = [None] * len(samfh.references)
    gene_trees = collections.defaultdict(lambda: IntervalTree())    
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
        # store gene by reference id in sam file
        gene_genome_map[gene_tid] = g
        # add gene to interval tree
        gene_interval = Interval(g.tx_start, g.tx_end, strand=g.strand, value=g.tx_name)
        gene_trees[chrom_tid].insert_interval(gene_interval)
    return gene_genome_map, gene_trees

def build_gene_to_genome_map(gene_feature_file, rname_prefix=None):
    # create arrays to map genes in bed file to genome 
    rname_prefix = '' if rname_prefix is None else rname_prefix
    bed_genome_map = {}    
    for g in GeneFeature.parse(open(gene_feature_file)):
        rname = rname_prefix + g.tx_name
        strand = 1 if g.strand == '-' else 0 
        exon_vectors = []
        if strand:
            exon_vectors = [(end, end-start) for start, end in g.exons]
            exon_vectors.reverse()
        else:
            exon_vectors = [(start, end-start) for start, end in g.exons]
        if rname in bed_genome_map:
            logging.error("Duplicate references %s found in bed file" % (rname))
        bed_genome_map[rname] = (g.chrom, strand, exon_vectors)
    return bed_genome_map

#def bed_to_genome_coords(line):
#    fields = line.split('\t')
#    chrom = fields[0]
#    tx_start = int(fields[1])
#    tx_end = int(fields[2])
#    name = fields[3]
#    strand = 1 if fields[5] == '-' else 0 
#    exon_count = int(fields[9])
#    exon_sizes = map(int, fields[10].split(',')[:-1])
#    exon_starts = map(int, fields[11].split(',')[:-1])
#    intervals = []
#    if strand:
#        exon_starts.reverse()
#        exon_sizes.reverse()
#        for e_start, e_size in itertools.izip(exon_starts, exon_sizes):
#            intervals.append((tx_start + e_start + e_size, e_size))
#    else:
#        for e_start, e_size in itertools.izip(exon_starts, exon_sizes):
#            intervals.append((tx_start + e_start, e_size))
#    return name, chrom, strand, intervals
#
#def parse_bed_genes(fileh):
#    for line in fileh:
#        if line is None:
#            continue
#        line = line.rstrip()
#        if line.startswith('#'):
#            logging.debug("BED skipping comment line: %s" % (line))
#            continue
#        if line.startswith('track'):
#            logging.debug("BED skipping track header line: %s"  % (line))
#            continue
#        yield bed_to_genome_coords(line)
#
#def build_gene_to_genome_map(bed_file, rname_prefix=None):
#    # create arrays to map genes in bed file to genome 
#    rname_prefix = '' if rname_prefix is None else rname_prefix
#    bed_genome_map = {}
#    for name,chrom,strand,intervals in parse_bed_genes(open(bed_file)):
#        rname = rname_prefix + name
#        if rname in bed_genome_map:
#            logging.error("Duplicate references %s found in bed file" % (rname))
#        bed_genome_map[rname] = (chrom, strand, intervals)
#    return bed_genome_map

def build_translation_table(samfh, gene_to_genome_map):

    new_header = samfh.header
    new_header['SQ'] = []
    rname_tid_map = {}
    current_tid = 0
    gene_table = []

    for tid, refdict in enumerate(samfh.header['SQ']):
        rname = refdict['SN']
        rlen = refdict['LN']
        # add lines to new header
        if rname not in gene_to_genome_map:
            new_header['SQ'].append(refdict)            
            rname_tid_map[rname] = current_tid
            current_tid += 1    
    # now that new header is built, can map references from 
    # the transcriptome mappings to genomic mappings
    for refdict in samfh.header['SQ']:
        rname = refdict['SN']
        rlen = refdict['LN']    
        # add lines to gene table and array
        if rname not in gene_to_genome_map:
            #logging.debug("Reference %s length=%d is a genomic reference with no transcriptome mapping" % (rname, rlen))
            gene_table.append((rname, rlen, rname_tid_map[rname], 0, None))
        else:
            chrom, strand, intervals = gene_to_genome_map[rname]
            if chrom not in rname_tid_map:
                logging.error("Translated rname=%s to chromosome=%s not found in translation table" % (rname, chrom))
                gene_table.append((rname, rlen, -1, 0, None))
            else:
                #logging.debug("Reference %s length=%d mapped to chromosome %s" % (rname, rlen, chrom))
                gene_table.append((rname, rlen, rname_tid_map[chrom], strand, intervals))
    return new_header, gene_table
    

def reverse_complement_MD_tag(val):
    x = 0
    mdops = []
    for y in xrange(len(val)):
        if val[y].isalpha():
            mdops.append(val[x:y])
            mdops.append(DNA_reverse_complement(val[y]))
            x = y + 1
    if x < len(val):
        mdops.append(val[x:])
    val = ''.join(mdops[::-1])
    return val

def translate_transcriptome_to_genomic_intervals(read, chrom, strand, intervals):
    # translate the read to genome space
    pos = 0
    read_pos = read.pos
    genomic_intervals = []
    for genomic_start,e_size in intervals:
        #logging.debug("read_pos=%d current_pos=%d exon=%d size=%d intervals=%s" % 
        #              (read.pos, pos, genomic_start, e_size, str(genomic_intervals)))
        # ending of read is before this exon
        if read.aend <= pos:
            break
        # beginning of read is in this exon
        if read_pos < pos + e_size:
            if strand:
                end = genomic_start - (read_pos - pos)
                start = genomic_start - min(e_size, read.aend - pos)
            else:
                start = genomic_start + read_pos - pos
                end = genomic_start + min(e_size, read.aend - pos)
            genomic_intervals.append((start, end))
            read_pos += (end - start)
        pos += e_size
    if strand:
        genomic_intervals.reverse()
    #logging.debug("read_pos=%d current_pos=%d exon=%d size=%d intervals=%s" % 
    #              (read.pos, pos, genomic_start, e_size, str(genomic_intervals)))
    return genomic_intervals

CIGAR_M = 0 #match  Alignment match (can be a sequence match or mismatch)
CIGAR_I = 1 #insertion  Insertion to the reference
CIGAR_D = 2 #deletion  Deletion from the reference
CIGAR_N = 3 #skip  Skipped region from the reference
CIGAR_S = 4 #softclip  Soft clip on the read (clipped sequence present in <seq>)
CIGAR_H = 5 #hardclip  Hard clip on the read (clipped sequence NOT present in <seq>)
CIGAR_P = 6 #padding  Padding (silent deletion from the padded reference sequence)
STRAND_FWD = 0
STRAND_REV = 1

def get_cigar(intervals):
    spliced = False
    cigar = [(CIGAR_M, intervals[0][1] - intervals[0][0])]
    prev_end = intervals[0][1]
    for start,end in intervals[1:]:
        cigar.append((CIGAR_N, start - prev_end))
        cigar.append((CIGAR_M, end - start))
        prev_end = end
        spliced = True
    #logging.debug("intervals: %s cigar: %s" % (str(intervals), str(cigar)))
    return spliced, cigar

def translate_read(read, chrom, strand, intervals):
    # skip unmapped reads
    if read.is_unmapped:
        return read
    elif chrom == -1:
        #logging.warning("discarded alignment %s that does not map to genomic references and cannot be translated" % (str(read)))
        # throw away reads that cannot be translated by
        # creating a dummy unmapped read
        a = pysam.AlignedRead()
        a.qname = read.qname
        a.seq = read.seq
        a.is_unmapped = True        
        a.rname = -1
        a.pos = -1
        a.mapq = 0
        a.mrnm = -1
        a.mpos = -1
        a.isize = 0
        a.qual = read.qual
        a.tags = [("XM", 0)]
        return a
    elif (chrom >= 0) and (intervals is None):
        # read maps directly to a genomic reference so simply
        # alter the reference id to correctly refer to the new 
        # SAM header
        read.rname = chrom
        return read   
    genomic_intervals = translate_transcriptome_to_genomic_intervals(read, chrom, strand, intervals)
    spliced, cigar = get_cigar(genomic_intervals)
    if spliced:
        read.tags = read.tags + [("XS", "-" if strand else "+")]
    # modify read
    read.rname = chrom
    read.pos = genomic_intervals[0][0]
    read.cigar = cigar
    # flip reads that aligned to negative strand genes        
    if strand == STRAND_REV:
        rev_quals = read.qual[::-1]
        read.is_reverse = not read.is_reverse
        read.seq = DNA_reverse_complement(read.seq)
        read.qual = rev_quals
        new_tags = []
        for name,val in read.tags:
            if name == 'MD':
                val = reverse_complement_MD_tag(val)
            new_tags.append((name,val))
        read.tags = new_tags
    return read
    
def translate_multihit_reads(insamfh, gene_table):
    # setup stats and debugging
    total_reads = 0
    total_alignments = 0
    unmapped = 0
    filtered = 0
    kept = 0
    debug_next, debug_step = 1e5,1e5
    # main loop to translate reads 
    for reads in parse_multihit_alignments(insamfh):
        total_reads += 1
        read_keys = set()
        unmapped_reads = []
        for read in reads:
            total_alignments += 1
            # read must have entry in gene lookup table
            rname_orig, rlen_orig, chrom, strand, intervals = gene_table[read.rname]
            #logging.debug("translating %s %s %d pos=%d chrom=%s strand=%s" % 
            #              (read.qname, rname_orig, rlen_orig, read.pos, chrom, strand))
            new_read = translate_read(read, chrom, strand, intervals)
            if new_read.is_unmapped:
                unmapped += 1
                unmapped_reads.append(new_read)
            else:
                k = (new_read.rname, new_read.pos, tuple(new_read.cigar), new_read.is_reverse)
                if k in read_keys:
                    filtered += 1
                    continue
                else:
                    read_keys.add(k)
                    kept += 1
                yield new_read
        if len(read_keys) == 0:
            yield unmapped_reads[0]
        if total_reads >= debug_next:
            logging.debug("Processed %d reads" % total_reads)
            logging.debug(" %d alignments" % total_alignments)
            logging.debug(" %d unmapped" % unmapped)
            logging.debug(" %d filtered transcriptome multihits" % filtered)
            logging.debug(" %d alignments used" % kept)
            debug_next += debug_step
    logging.info("Processed %d reads" % total_reads)
    logging.info(" %d alignments" % total_alignments)
    logging.info(" %d unmapped" % unmapped)
    logging.info(" %d filtered transcriptome multihits" % filtered)
    logging.info(" %d alignments used" % kept)

def transcriptome_to_genome(input_sam_file, output_sam_file, gene_to_genome_map): 
    insamfh = pysam.Samfile(input_sam_file, "r")
    new_header, gene_table = build_translation_table(insamfh, gene_to_genome_map)
    outsamfh = pysam.Samfile(output_sam_file, "wh", header=new_header)
    for read in translate_multihit_reads(insamfh, gene_table):
        outsamfh.write(read)
    outsamfh.close()
    insamfh.close()
    #if strand and len(cigar) > 2:
    #    print read.seq
    #    print insamfh.getrname(read.rname), read.pos, read.aend, read.cigar
    #    #print 'old', insamfh.getrname(read.rname), read.pos, read.cigar
    #    #print 'new', insamfh.getrname(chrom), genome_pos, cigar 


def main():
    # TODO: switch from argparse to OptionParser for backwards compatibility
    import argparse
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")    
    parser = argparse.ArgumentParser()
    parser.add_argument('--build', action="store_true", default=False)
    parser.add_argument('--rname-prefix', dest="rname_prefix", default=None)
    parser.add_argument('--genes', dest="gene_feature_file", default=None)
    parser.add_argument('map_file')
    parser.add_argument('sam_file', nargs="?")
    options = parser.parse_args()

    if options.build:
        if options.gene_feature_file is None:
            parser.error("--build requires gene feature file (--genes=file)")
        gene_to_genome_map = build_gene_to_genome_map(options.gene_feature_file, options.rname_prefix)
        pickle.dump(gene_to_genome_map, open(options.map_file, "w"))
    else:
        gene_to_genome_map = pickle.load(open(options.map_file))
        transcriptome_to_genome(options.sam_file, "-", gene_to_genome_map)

if __name__ == '__main__': main()
