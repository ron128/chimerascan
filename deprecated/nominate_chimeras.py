'''
Created on Jan 11, 2011

@author: mkiyer
'''
'''
Created on Jan 5, 2011

@author: mkiyer
'''
import collections
import logging
import operator
import os

import pysam

from bx.intersection import Interval, IntervalTree
from base import parse_multihit_alignments, get_aligned_read_intervals
from feature import GeneFeature

def get_unpaired_reads(bamfh):
    pe_reads = ([], [])
    # reads must be binned by qname, mate, hit, and segment
    # so initialize to mate 0, hit 0, segment 0
    num_reads = 0
    prev_qname = None
    for read in bamfh:
        # if read is mapped in proper pair continue
        if read.is_proper_pair:
            continue
        # get read attributes
        qname = read.qname
        mate = 0 if read.is_read1 else 1
        hit_ind = read.opt('HI')
        num_hits = read.opt('IH')
        seg_ind = read.opt('XI')
        num_segs = read.opt('XN')
        # if query name changes we have completely finished
        # the fragment and can reset the read data
        if num_reads > 0 and qname != prev_qname:
            yield pe_reads
            # reset state variables
            pe_reads = ([], [])
            num_reads = 0
        prev_qname = qname
        # initialize mate hits
        if len(pe_reads[mate]) == 0:
            pe_reads[mate].extend([list() for x in xrange(num_hits)])
        mate_reads = pe_reads[mate]
        # initialize hit segments
        if len(mate_reads[hit_ind]) == 0:
            mate_reads[hit_ind].extend([None for x in xrange(num_segs)])
        # add segment to hit/mate/read
        mate_reads[hit_ind][seg_ind] = read
        num_reads += 1
    if num_reads > 0:
        yield pe_reads

def map_reads_to_references(pe_reads, tid_list):
    # bin reads by reference name to find reads that pairs
    # to the same gene/chromosome
    ref_dict = collections.defaultdict(lambda: [[], []])
    partial_mapping_dict = collections.defaultdict(lambda: [[], []])
    nonmapping_dict = collections.defaultdict(lambda: [[], []])
    for mate, mate_hits in enumerate(pe_reads):
        for hitsegs in mate_hits:
            num_unmapped_segs = sum(1 if r.is_unmapped else 0 for r in hitsegs)
            if num_unmapped_segs == len(hitsegs):
                assert len(hitsegs) == 1
                nonmapping_dict[mate].append(hitsegs[0])
            elif num_unmapped_segs > 0:
                partial_mapping_dict[mate].append(hitsegs)
            else:
                assert len(hitsegs) == 1
                read = hitsegs[0]
                mate_pairs = ref_dict[read.rname]
                mate_pairs[mate].append(read)
    return ref_dict, partial_mapping_dict, nonmapping_dict




def find_discordant_mappings():
    # if there are any unmapped segments, this read cannot be considered
    # concordant, but it might be a putative chimeric read where
    # the chimeric junction creates an unmapped segment
    if len(unmapped_segs) > 0:
        #yield make_unmapped_read(mate, 
        return

def parse_pe_multihit_alignments(samfh, 
                                 remove_unmapped=False, 
                                 contam_tids=None):
    if contam_tids is None:
        contam_tids = set()
    for read_hits in parse_multihit_alignments(samfh):
        mate_reads = ([], [])
        for read in read_hits:
            if remove_unmapped and read.is_unmapped:
                continue 
            if read.rname in contam_tids:
                continue            
            mate = 0 if read.is_read1 else 1
            mate_reads[mate].append(read)
        yield mate_reads

def select_best_reads(reads):
    # sort reads by number of mismatches
    mismatch_read_tuples = sorted([(r.opt('NM'), r) for r in reads], key=operator.itemgetter(0))
    # set this to allow 'best reads' to have unequal 
    # amounts of mismatches
    mismatch_tolerance = 0
    best_mismatches = mismatch_read_tuples[0][0]
    best_reads = [mismatch_read_tuples[0][1]]
    # choose reads within a certain mismatch tolerance
    for mismatches, read in mismatch_read_tuples:
        if mismatches > best_mismatches + mismatch_tolerance:
            break
        best_reads.append(read)
    return best_reads
#    # weight best reads by expression level
#    expr_read_tuples = []
#    for read in best_reads:
#        max_expr = 0
#        for start, end in get_aligned_read_intervals(reads):
#            # ensure exon entirely encompasses the read
#            exon_ids = [hit.value for hit in exon_trees[read.rname].find(start, end)
#                        if hit.start <= start and hit.end >= end]
#            # extract expression levels
#            for exon_id in exon_ids:
#                exon_info_tuples = exon_id_map[exon_id]
#                for exon_info_tuple in exon_info_tuples:
#                    if exon_info_tuple[2] > max_expr:
#                        max_expr = exon_info_tuple[2]
#        expr_read_tuples.append((max_expr, read))

def discordant_reads_to_chimeras(samfh, contam_tids):    
    # establish counters
    total_reads = 0
    total_alignment_hits = 0
    both_non_mapping = 0
    single_non_mapping = 0
    split_mates = 0    
    # logging output
    debug_every = 1e5
    debug_next = debug_every    
    for mate_hits in parse_pe_multihit_alignments(samfh, 
                                                  remove_unmapped=True, 
                                                  contam_tids=contam_tids):
        # logging debug output
        total_reads += 1
        total_alignment_hits += len(mate_hits[0]) + len(mate_hits[1])
        if total_reads == debug_next:
            debug_next += debug_every
            logging.debug("Processed reads=%d alignments=%d" % 
                          (total_reads, total_alignment_hits))
        read1_hits, read2_hits = mate_hits
        if len(read1_hits) == 0 and len(read2_hits) == 0:
            both_non_mapping += 1
        elif len(read1_hits) == 0 or len(read2_hits) == 0:
            single_non_mapping += 1
        else:
            split_mates += 1
            read1_best_hits = select_best_reads(read1_hits)
            read2_best_hits = select_best_reads(read2_hits)
            yield read1_best_hits, read2_best_hits
    logging.info("Total reads=%d" % (total_reads))
    logging.info("Total alignment hits=%d" % (total_alignment_hits))
    logging.info("Both non-mapping pairs=%d" % (both_non_mapping))
    logging.info("Single non-mapping pairs=%d" % (single_non_mapping))
    logging.info("Chimeric reads=%d" % (split_mates))

def build_exon_expression_map(expression_file):
    # build gene expression lookup
    exon_expr_map = {}
    for line in open(expression_file):
        fields = line.strip().split('\t')
        if line.startswith('#'):
            continue
        gene_name = fields[0]
        if gene_name not in exon_expr_map:
            exon_expr_map[gene_name] = {}
        exon_num = int(fields[1])        
        expr = float(fields[7]) + float(fields[8])
        exon_expr_map[gene_name][exon_num] = expr
    for gene_name in exon_expr_map.iterkeys():
        sorted_exon_exprs = [v for k,v in sorted(exon_expr_map[gene_name].iteritems(), key=operator.itemgetter(0))]        
        exon_expr_map[gene_name] = sorted_exon_exprs
    return exon_expr_map
    
def build_genome_gene_map(samfh, genefile):
    rname_tid_map = dict((rname,i) for i,rname in enumerate(samfh.references))    
    trees = collections.defaultdict(lambda: IntervalTree())    
    exon_key_map = {}
    exon_id_map = collections.defaultdict(lambda: [])
    current_exon_id = 0
    # build genome -> gene lookup
    for g in GeneFeature.parse(open(genefile)):
        if g.chrom not in rname_tid_map:
            continue        
        tid = rname_tid_map[g.chrom]
        for exon_num, exon in enumerate(g.exons):
            start, end = exon
            exon_key = (g.chrom, start, end, g.strand)
            if exon_key in exon_key_map:
                exon_id = exon_key_map[exon_key]
            else:
                exon_id = current_exon_id
                current_exon_id += 1
                exon_key_map[exon_key] = exon_id
                trees[tid].insert_interval(Interval(start, end, 
                                                    strand=g.strand, chrom=g.chrom, value=exon_id))
            exon_id_map[exon_id].append((g, exon_num))
            #print 'exon id', exon_id, 'genes', [(x[0].name,x[0].genomic_exons[x[1]]) for x in intervals[exon_id]]
    return trees, exon_id_map

def find_read_exons(read, exon_trees, exon_id_map):
    # get the distinct genomic intervals that the
    # read aligns to
    intervals = get_aligned_read_intervals(read)
    total_length = sum((end-start) for start,end in intervals)     
    read_exon_data = []
    start, end = intervals[0]
    exon_ids = set(hit.value for hit in exon_trees[read.rname].find(start, end)
                   if start >= hit.start and end <= hit.end)
    for start,end in intervals[1:]:
        # find exons that contain this read interval
        exon_ids.intersection_update(hit.value for hit in exon_trees[read.rname].find(start, end)
                                     if start >= hit.start and end <= hit.end)

        for hit in exon_trees[read.rname].find(start, end):
            if start >= hit.start and end <= hit.end:
                g, exon_num = exon_id_map[hit.value]        
                read_exon_data.append((start, end, weight, g, exon_num))
    return read_exon_data

def read_pair_to_gene_bedpe(read1, read2, exon_trees, exon_id_map):
    read1_exon_data = find_read_exons(read1)
    read2_exon_data = find_read_exons(read2)
    for r1_start, r1_end, r1_weight, r1_gene, r1_exon_num in read1_exon_data:
        for r2_start, r2_end, r2_weight, r2_gene, r2_exon_num in read2_exon_data:
            score = float(r1_weight + r2_weight) / 2.0

    #chr8    74203462        74203486        chr17   1326672 1326696 PATHBIO-SOLEXA2_30TUEAAXX:3:1:962:931   1       -1      -1      MULTI;CTCCATGCAGATGATGCCGTATTTA;AAAAAATATCACAGGTAACAGAATG       chr8    74203277        74206379        uc003xzh.1      0       -       74203278        74204926        0       6       210,110,138,167,109,279,        0,509,730,1196,1646,2823,


def find_discordant_gene_pairs(pe_reads, tid_list, min_isize, max_isize,
                               library_type):    
    same_strand = (library_type[0] == library_type[1])
    # check for mapping to same gene within insert size range
    gene_pairs = []
    genome_pairs = []
    for rname, mate_pairs in ref_dict.iteritems():
        pairs_list = genome_pairs if tid_list[rname] == None else gene_pairs
        if len(mate_pairs[0]) > 0 and len(mate_pairs[1]) > 0:
            # ensure distance is within insert size range
            # and strandedness matches library type
            for r1 in mate_pairs[0]:
                for r2 in mate_pairs[1]:
                    # check strandedness of pairs
                    if same_strand != (r1.is_reverse == r2.is_reverse):
                        continue
                    # check insert size                                         
                    if r1.pos > r2.pos:
                        isize = r1.aend - r2.pos
                    else:
                        isize = r2.aend - r1.pos
                    if isize < min_isize or isize > max_isize:
                        continue
                    strand_matches_library = (r1.is_reverse != library_type[0] or
                                              r2.is_reverse != library_type[1])
                    # this is a concordant read pair
                    pairs_list.append((r1, r2, strand_matches_library))
    return gene_pairs, genome_pairs


def get_tids(samfh, rnames):
    rname_tid_map = dict((rname,i) for i,rname in enumerate(samfh.references))    
    return [rname_tid_map[rname] for rname in rnames]

def nominate_chimeras(bam_file, gene_file, expression_file, 
                      library_type='fr', 
                      contam_refs=None):

    logging.info("Reading chimerascan gene index")
    bamfh = pysam.Samfile(input_bam_file, "rb")
    gene_file = os.path.join(options.index_dir, config.GENE_FEATURE_FILE)
    gene_genome_map = build_gene_genome_map(bamfh, gene_file)

    
    bamfh = pysam.Samfile(bam_file, "rb")
    # get contaminant reference tids
    if contam_refs is None:
        contam_tids = set()
    else:
        contam_tids = set(get_tids(bamfh, contam_refs))
    # build interval trees from annotated exons
    logging.debug("Building gene lookup tables...")
    exon_trees, exon_intervals = build_genome_gene_map(bamfh, gene_file, expression_file)
    # load expression data
    logging.debug("Reading gene expression data...")
    exon_expr_map = build_exon_expression_map(expression_file)
    # filter out contaminant reads and prioritize 
    # mates by number of mismatches
    for read1_hits, read2_hits in discordant_reads_to_chimeras(bamfh, contam_tids):    
        read1_num_hits = len(read1_hits)
        read2_num_hits = len(read2_hits)
        multimap = (read1_num_hits > 1) or (read2_num_hits > 1)                
        for r1 in read1_hits:
            for r2 in read2_hits:
                yield r1,r2,multimap

    bamfh.close()

    #
    # Step 2: Compare BEDPE using BEDTools to all genes
    #
    bedpe_overlap_file = os.path.join(tmp_dir, "chimera_gene_overlap_bedpe.txt")
    if not os.path.exists(bedpe_overlap_file):
        logging.info("%s: Finding overlapping genes" % (job_name))
        args = [pairtobed_bin, "-type", "both", 
                "-a", prev_output_file, "-b", gene_bed_file]
        logging.debug("%s: args=%s" % (job_name, ' '.join(args)))
        f = open(bedpe_overlap_file, "w")
        if subprocess.call(args, stdout=f) != JOB_SUCCESS:
            logging.error("%s: Error finding overlapping genes" % (job_name))    
            return JOB_ERROR
        f.close()
    prev_output_file = bedpe_overlap_file

def build_gene_genome_map(samfh, genefile):
    rname_tid_map = dict((rname,i) for i,rname in enumerate(samfh.references))
    gene_genome_map = [None] * len(samfh.references)
    # build gene -> genome lookup
    for g in GeneFeature.parse(open(genefile)):
        name = config.GENE_REF_PREFIX + g.tx_name
        if name not in rname_tid_map:
            continue
        if g.chrom not in rname_tid_map:
            continue
        gene_tid = rname_tid_map[name]
        chrom_tid = rname_tid_map[g.chrom]        
        gene_genome_map[gene_tid] = (chrom_tid, g.tx_start, g.tx_end, g.strand)
    return gene_genome_map


    #logging.debug("Clustering chimeric gene pairs")
    #exon_pair_reads = cluster_reads(bamfh, exon_trees, exon_intervals, contam_tids, options.min_dist)    
    #logging.debug("Filtering chimeric genes")
    #filter_chimeric_genes(exon_pair_reads, exon_intervals)  

    
def main():
    from optparse import OptionParser
    import sys
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = OptionParser("usage: %prog [options] <bam> <out.bedpe>")
    parser.add_option('--max-fragment-length', dest="max_fragment_length", 
                      type="int", default=1000)
    parser.add_option('--library-type', dest="library_type", default="fr")
    parser.add_option("--index", dest="index_dir",
                      help="Path to chimerascan index directory")    
    options, args = parser.parse_args()
    input_bam_file = args[0]
    output_bam_file = args[1]
    
    gene_file = os.path.join(config)
    
    library_type = parse_library_type(options.library_type)
    logging.info("Finding discordant reads")
    logging.debug("Input file: %s" % (input_bam_file))
    logging.debug("Output file: %s" % (output_bam_file))
    logging.debug("Library type: %s" % (library_type))


    logging.info("Checking reads")
    outfh = pysam.Samfile(output_bam_file, "wb")
    find_discordant_reads(bamfh, outfh, 
                          options.min_fragment_length,
                          options.max_fragment_length,
                          library_type)
    
if __name__ == '__main__':
    main()