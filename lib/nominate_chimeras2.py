'''
Created on Jan 11, 2011

@author: mkiyer
'''
import collections
import logging
import operator
import os
import sys

import pysam

from bx.intersection import Interval, IntervalTree
from feature import GeneFeature
import config
from base import parse_library_type

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

def map_reads_to_references(pe_reads):
    # bin reads by reference name to find reads that pairs
    # to the same gene/chromosome
    ref_dict = collections.defaultdict(lambda: [[], []])
    for mate, mate_hits in enumerate(pe_reads):
        for hitsegs in mate_hits:
            if any(r.is_unmapped for r in hitsegs):
                continue
            assert len(hitsegs) == 1
            read = hitsegs[0]
            mate_pairs = ref_dict[read.rname]
            mate_pairs[mate].append(read)
    return ref_dict


def interval_overlap(chrom1, start1, end1, chrom2, start2, end2):
    if chrom1 != chrom2:
        return False
    return (start1 < end2) and (start2 < end1)    

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
    chrom1, start1, end1, strand1 = fiveprime_gene.chrom, fiveprime_gene.tx_start, fiveprime_gene.tx_end, fiveprime_gene.strand
    chrom2, start2, end2, strand2 = threeprime_gene.chrom, threeprime_gene.tx_start, threeprime_gene.tx_end, threeprime_gene.strand
    # interchromosomal
    if chrom1 != chrom2:
        return CHIMERA_INTERCHROMOSOMAL, None    
    # orientation
    same_strand = strand1 == strand2
    # genes on same chromosome so check overlap
    is_overlapping = (start1 < end2) and (start2 < end1)            
    if is_overlapping:
        if not same_strand:
            if ((start1 <= start2 and strand1 == "+") or
                (start1 > start2 and strand1 == "-")):                    
                return (CHIMERA_OVERLAP_CONVERGE, 0)
            elif ((start1 >= start2 and strand1 == "+") or
                  (start1 < start2 and strand1 == "-")):
                return (CHIMERA_OVERLAP_DIVERGE, 0)
        else:
            if ((start1 <= start2 and strand1 == "+") or
                (start1 > start2 and strand1 == "-")):
                return (CHIMERA_OVERLAP_SAME, 0)
            elif ((start1 >= start2 and strand1 == "+") or
                  (start1 < start2 and strand1 == '-')):
                return (CHIMERA_OVERLAP_COMPLEX, 0)    
    # if code gets here then the genes are on the same chromosome but do not
    # overlap.  first calculate distance (minimum distance between genes)
    if start1 <= start2:
        distance = start2 - end1
        between_interval = (end1, start2)
    else:
        distance = start1 - end2
        between_interval = (end2, start1)
    # check whether there are genes intervening between the
    # chimera candidates
    genes_between = []
    genes_between_same_strand = []
    for hit in gene_trees[chrom1].find(*between_interval):
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
            if ((start1 <= start2 and strand1 == "+") or
                (start1 > start2 and strand1 == "-")):                    
                return (CHIMERA_ADJ_CONVERGE, distance)
            elif ((start1 >= start2 and strand1 == "+") or
                  (start1 < start2 and strand1 == "-")):
                return (CHIMERA_ADJ_DIVERGE, distance)
            elif ((start1 <= start2 and strand1 == "+") or
                  (start1 > start2 and strand1 == "-")):
                return (CHIMERA_ADJ_SAME, distance)
            elif ((start1 >= start2 and strand1 == "+") or
                  (start1 < start2 and strand1 == '-')):
                return (CHIMERA_ADJ_COMPLEX, distance)
        else:
            return CHIMERA_INTRA_COMPLEX, distance    
    return CHIMERA_UNKNOWN, distance


def select_best_segments(seg_list):
    # sort by number of unmapped segments
    sorted_segs = []
    for segs in seg_list:
        num_unmapped_segs = sum(1 if r.is_unmapped else 0 for r in segs)
        sorted_segs.append((num_unmapped_segs, segs))
    sorted_segs = sorted(sorted_segs, key=operator.itemgetter(0))
    best_num_unmapped_segs = sorted_segs[0][0]
    best_segs = []
    for num_unmapped_segs, segs in sorted_segs:
        if num_unmapped_segs > best_num_unmapped_segs:
            break
        best_segs.append(segs)
    return best_segs

def select_best_mismatches(seg_list, mismatch_tolerance=0):
    if len(seg_list) == 0:
        return []
    sorted_segs_by_mismatch = []
    for segs in seg_list:
        # keep track of number of alignment mismatches in the pair
        num_mismatches = sum(r.opt('NM') for r in segs
                             if not r.is_unmapped)
        sorted_segs_by_mismatch.append((num_mismatches, segs))
    sorted_segs_by_mismatch = sorted(sorted_segs_by_mismatch, key=operator.itemgetter(0))    
    best_num_mismatches = sorted_segs_by_mismatch[0][0]
    # choose reads within a certain mismatch tolerance
    best_segs = []
    for mismatches, segs in sorted_segs_by_mismatch:
        if mismatches > best_num_mismatches + mismatch_tolerance:
            break
        best_segs.append(segs)
    return best_segs


def find_discordant_pairs(pe_reads, tid_list, genome_tid_set, 
                          gene_genome_map, gene_trees, max_isize, 
                          contam_tids):    
    unmapped_hits = ([],[])
    genomic_hits = ([],[])
    mapped_hits = ([],[])
    # classify reads and mapping/nonmapping and determine
    # set of reference names for read1 and read2
    for mate, mate_hits in enumerate(pe_reads):
        for segs in mate_hits:
            if len(segs) == 1 and segs[0].is_unmapped:
                # bin nonmapping segments separately
                unmapped_hits[mate].append(segs[0])
            elif any(tid_list[seg.rname] is None for seg in segs):
                # ignore genomic mappings
                genomic_hits[mate].append(segs)
            else:
                # keep mapped segments
                mapped_hits[mate].append(segs)
    
    # filter hits based on library type and bin by
    # 5' and 3' based on strand matching library type
    filtered_hits = (([], []), ([], []))
    for mate, mate_hits in enumerate(mapped_hits):
        for segs in mate_hits:
            #print 'num_unmapped', num_unmapped_segs, 'mate', mate, '|||'.join(['%s:%d' % (tid_list[r.rname], r.pos) for r in segs])
            # determine sense/antisense (5' or 3')
            sense_index = 0
            for seg in segs:
                # ignore unmapped segments
                if seg.is_unmapped:
                    continue
                # TODO: support multiple library types here
                sense_index = 1 if seg.is_reverse else 0
                break
            filtered_hits[sense_index][mate].append(segs)

    # join 5' and 3' hits
    concordant_partially_mapped = []
    discordant_unmapped = []
    discordant_mapped = []    
    for mate1, mate2 in ((0, 1), (1, 0)):
        for fiveprime_segs in filtered_hits[0][mate1]:
            for threeprime_segs in filtered_hits[1][mate2]:
                # combine 5'/3' partners
                segs = []
                segs.extend(fiveprime_segs)
                segs.extend(threeprime_segs)
                # check rnames
                rnames = set(seg.rname for seg in segs
                             if not seg.is_unmapped)
                assert len(rnames) > 0
                if len(rnames) == 1:
                    # only one gene reference so this is a 
                    # concordant pair
                    concordant_partially_mapped.append(segs)
                elif len(rnames) > 2:
                    discordant_unmapped.append(segs)
                else:
                    # double check that the first/last segments are mapped
                    if segs[0].is_unmapped or segs[-1].is_unmapped:                    
                        discordant_unmapped.append(segs)
                    else:
                        discordant_mapped.append(segs)

    # select hits with fewest mismatches, keeping categories separate for now
    concordant_partially_mapped = select_best_mismatches(concordant_partially_mapped)
    discordant_unmapped = select_best_mismatches(discordant_unmapped)
    discordant_mapped = select_best_mismatches(discordant_mapped)
    
    #return discordant_mapped, discordant_unmapped

    for segs in discordant_mapped:
        print '|||'.join(['%s:%d' % ('None' if r.rname == -1 else tid_list[r.rname], r.pos) for r in segs])

    #print 'mate', mate1, '|||'.join(['%s:%d' % (tid_list[r.rname], r.pos) for r in fiveprime_segs])
    #print 'mate', mate2, '|||'.join(['%s:%d' % (tid_list[r.rname], r.pos) for r in threeprime_segs])
    #print fiveprime_segs
    #print threeprime_segs
               



def get_gene_tids(bamfh):
    gene_tids = []
    for ref in bamfh.references:
        if ref.startswith(config.GENE_REF_PREFIX):
            gene_tids.append(ref)
        else:
            gene_tids.append(None)
    return gene_tids

def get_genome_tids(bamfh):
    tids = set()
    for tid,ref in enumerate(bamfh.references):
        if not ref.startswith(config.GENE_REF_PREFIX):
            tids.add(tid)
    return tids

def get_tids(samfh, rnames):
    rname_tid_map = dict((rname,i) for i,rname in enumerate(samfh.references))    
    contam_tids = []
    for rname in rnames:
        if rname not in rname_tid_map:
            logging.warning("Reference %s not found in SAM file.. ignoring" % (rname))
        else:
            contam_tids.append(rname_tid_map[rname])
    return contam_tids

def nominate_chimeras(input_bam_file, output_bedpe_file, gene_file,
                      max_isize, library_type='fr', contam_refs=None):
    logging.info("Nominating chimeras")
    logging.debug("Input file: %s" % (input_bam_file))
    logging.debug("Output file: %s" % (output_bedpe_file))
    logging.debug("Library type: %s" % (library_type))
    logging.debug("Contaminant references: %s" % (contam_refs))

    logging.info("Read BAM file")
    bamfh = pysam.Samfile(input_bam_file, "rb")    
    # build a map of gene name to genome coords
    logging.info("Reading gene index")
    gene_genome_map, gene_trees = build_gene_maps(bamfh, gene_file)
    # parse the library type string into a tuple
    library_type = parse_library_type(library_type)
    assert library_type == (0, 1)    
    #same_strand = (library_type[0] == library_type[1])    
    # get contaminant reference tids
    if contam_refs is None:
        contam_tids = set()
    else:
        contam_tids = set(get_tids(bamfh, contam_refs))
    # search for discordant reads that represent chimeras
    logging.info("Finding discordant pairs")    
    gene_tid_list = get_gene_tids(bamfh)
    genome_tid_set = get_genome_tids(bamfh)
    for pe_reads in get_unpaired_reads(bamfh):        
        find_discordant_pairs(pe_reads, gene_tid_list, genome_tid_set, gene_genome_map, 
                              gene_trees, max_isize, contam_tids)

def main():
    from optparse import OptionParser
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = OptionParser("usage: %prog [options] <bam> <out.bedpe>")
    parser.add_option('--max-fragment-length', dest="max_fragment_length", 
                      type="int", default=1000)
    parser.add_option('--library-type', dest="library_type", default="fr")
    parser.add_option("--index", dest="index_dir",
                      help="Path to chimerascan index directory")
    parser.add_option("--contam-refs", dest="contam_refs", default=None)
    options, args = parser.parse_args()
    input_bam_file = args[0]
    output_bedpe_file = args[1]
    gene_file = os.path.join(options.index_dir, config.GENE_FEATURE_FILE)
    nominate_chimeras(input_bam_file, output_bedpe_file, gene_file, 
                      options.max_fragment_length, options.library_type,
                      options.contam_refs)

if __name__ == '__main__':
    main()
    
#def map_reads_to_references(pe_reads, tid_list):
#    # bin reads by reference name to find reads that pairs
#    # to the same gene/chromosome
#    ref_dict = collections.defaultdict(lambda: [[], []])
#    partial_mapping_dict = collections.defaultdict(lambda: [[], []])
#    nonmapping_dict = collections.defaultdict(lambda: [[], []])
#    for mate, mate_hits in enumerate(pe_reads):
#        for hitsegs in mate_hits:
#            num_unmapped_segs = sum(1 if r.is_unmapped else 0 for r in hitsegs)
#            if num_unmapped_segs == len(hitsegs):
#                assert len(hitsegs) == 1
#                nonmapping_dict[mate].append(hitsegs[0])
#            elif num_unmapped_segs > 0:
#                partial_mapping_dict[mate].append(hitsegs)
#            else:
#                assert len(hitsegs) == 1
#                read = hitsegs[0]
#                mate_pairs = ref_dict[read.rname]
#                mate_pairs[mate].append(read)
#    return ref_dict, partial_mapping_dict, nonmapping_dict

#def discordant_reads_to_chimeras(samfh, contam_tids):    
#    # establish counters
#    total_reads = 0
#    total_alignment_hits = 0
#    both_non_mapping = 0
#    single_non_mapping = 0
#    split_mates = 0    
#    # logging output
#    debug_every = 1e5
#    debug_next = debug_every    
#    for mate_hits in parse_pe_multihit_alignments(samfh, 
#                                                  remove_unmapped=True, 
#                                                  contam_tids=contam_tids):
#        # logging debug output
#        total_reads += 1
#        total_alignment_hits += len(mate_hits[0]) + len(mate_hits[1])
#        if total_reads == debug_next:
#            debug_next += debug_every
#            logging.debug("Processed reads=%d alignments=%d" % 
#                          (total_reads, total_alignment_hits))
#        read1_hits, read2_hits = mate_hits
#        if len(read1_hits) == 0 and len(read2_hits) == 0:
#            both_non_mapping += 1
#        elif len(read1_hits) == 0 or len(read2_hits) == 0:
#            single_non_mapping += 1
#        else:
#            split_mates += 1
#            read1_best_hits = select_best_reads(read1_hits)
#            read2_best_hits = select_best_reads(read2_hits)
#            yield read1_best_hits, read2_best_hits
#    logging.info("Total reads=%d" % (total_reads))
#    logging.info("Total alignment hits=%d" % (total_alignment_hits))
#    logging.info("Both non-mapping pairs=%d" % (both_non_mapping))
#    logging.info("Single non-mapping pairs=%d" % (single_non_mapping))
#    logging.info("Chimeric reads=%d" % (split_mates))
