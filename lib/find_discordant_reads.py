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
from base import parse_library_type
from seq import DNA_reverse_complement
from gene_to_genome import build_gene_maps, get_gene_tids
from alignment_parser import parse_unpaired_reads

# constants
ChimeraMate = collections.namedtuple('ChimeraMate',
                                     ['interval', 'seq', 'qual', 'is_spanning'])

# mapping codes
NM = 0
MULTIMAP = 1
GENOME = 2
MAP = 3
    
class Chimera(object):
    DISCORDANT_INNER = 0
    DISCORDANT_5P = 1
    DISCORDANT_3P = 2
    DISCORDANT_OVERLAPPING = 3
    NONMAPPING = 4
    _discordant_types = ["DISCORDANT_INNER",
                         "DISCORDANT_5PRIME",
                         "DISCORDANT_3PRIME",
                         "DISCORDANT_OVERLAPPING",
                         "NONMAPPING"]
    
    def __init__(self, qname, discordant_type, mate5p, mate3p, 
                 read1_is_5prime):
        self.qname = qname
        self.discordant_type = discordant_type
        self.mate5p = mate5p
        self.mate3p = mate3p
        self.read1_is_5prime = read1_is_5prime
    
    def to_bedpe(self):
        s = [self.mate5p.interval.chrom,
             self.mate5p.interval.start,
             self.mate5p.interval.end,
             self.mate3p.interval.chrom,
             self.mate3p.interval.start,
             self.mate3p.interval.end,
             self.qname,
             1,
             self.mate5p.interval.strand,
             self.mate3p.interval.strand,
             self.mate5p.seq,
             self.mate3p.seq,
             self.mate5p.qual,
             self.mate3p.qual,
             1 if self.mate5p.is_spanning else 0,
             1 if self.mate3p.is_spanning else 0,
             Chimera._discordant_types[self.discordant_type],
             1 if self.read1_is_5prime else 0]             
        return '\t'.join(map(str, s))

    @staticmethod
    def from_bedpe(line):
        fields = line.strip().split('\t')
        interval5p = Interval(int(fields[1]),
                              int(fields[2]),
                              chrom=fields[0],
                              strand=fields[8])
        mate5p = ChimeraMate(interval5p, fields[10], fields[12], bool(int(fields[14])))
        interval3p = Interval(int(fields[4]),
                              int(fields[5]),
                              chrom=fields[3],
                              strand=fields[9])
        mate3p = ChimeraMate(interval3p, fields[11], fields[13], bool(int(fields[15])))
        discordant_type = Chimera._discordant_types.index(fields[16])        
        read1_is_5prime = bool(int(fields[17]))        
        return Chimera(fields[6], discordant_type, mate5p, mate3p, read1_is_5prime)



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


def select_best_mismatches(reads, mismatch_tolerance=0):
    if len(reads) == 0:
        return []
    # sort reads by number of mismatches
    mapped_reads = []
    unmapped_reads = []
    for r in reads:
        if r.is_unmapped:
            unmapped_reads.append(r)
        else:
            mapped_reads.append((r.opt('NM'), r))
    if len(mapped_reads) == 0:
        return unmapped_reads
    sorted_reads = sorted(mapped_reads, key=operator.itemgetter(0))
    best_nm = sorted_reads[0][0]
    worst_nm = sorted_reads[-1][0]
    sorted_reads.extend((worst_nm, r) for r in unmapped_reads)
    # choose reads within a certain mismatch tolerance
    best_reads = []
    for mismatches, r in sorted_reads:
        if mismatches > best_nm + mismatch_tolerance:
            break
        best_reads.append(r)
    return best_reads

def process_partition(partitions):
    for partition in partitions:
        splits5p = []
        splits3p = []
        codes5p = set()
        codes3p = set()
        for split_reads in partition:
            reads5p = []
            reads3p = []
            # TODO: we select reads in the best 'strata', that is, the
            # set of reads with fewest mismatches to the reference.  is
            # this the best strategy, or should other metrics be employed
            # to choose from among multimapping reads?
            best_reads = select_best_mismatches(split_reads)
            for r in best_reads:
                if r.is_unmapped:
                    reads3p.append(r)
                    reads5p.append(r)
                    if r.rname != -1:
                        code = GENOME
                    elif r.opt('XM') > 0:
                        code = MULTIMAP
                    else:
                        code = NM
                    codes5p.add(code)
                    codes3p.add(code)
                else:
                    if r.is_reverse:
                        # determine sense/antisense by assuming that
                        # 5' reads are sense and 3' reads are antisense                    
                        reads3p.append(r)
                        codes3p.add(MAP)
                    else:
                        reads5p.append(r)
                        codes5p.add(MAP)
            splits5p.append(reads5p)
            splits3p.append(reads3p)
        if all(len(reads) > 0 for reads in splits5p):
            yield 0, codes5p, splits5p
        if all(len(reads) > 0 for reads in splits3p):
            yield 1, codes3p, splits3p

def find_first_split(splits):
    rname_dict = collections.defaultdict(lambda: collections.defaultdict(lambda: []))
    split_ind = 0
    mapped_inds = []
    # find first split with mapped reads and 
    # save reads associated with reference names
    for ind,split_reads in enumerate(splits):
        split_rname_dict = collections.defaultdict(lambda: [])
        for r in split_reads:
            if r.is_unmapped:
                continue
            split_rname_dict[r.rname].append(r)
        if len(split_rname_dict) > 0:
            if len(rname_dict) == 0:   
                rnames = set(split_rname_dict)
            else:
                rnames = set(rname_dict).intersection(split_rname_dict)
                if len(rnames) == 0:
                    # if no ref names in common found a split so return
                    break
                # remove all ref names that are not in the intersection
                remove_rnames = set(rname_dict).difference(rnames)
                for rname in remove_rnames:
                    del rname_dict[rname]
            # add reads from this split
            for rname in rnames:
                # TODO: what happens with >1 mapping from split segment 
                # to same gene?
                rname_dict[rname][ind].extend(split_rname_dict[rname])
            mapped_inds.append(ind)            
        split_ind += 1
    return rname_dict, mapped_inds, split_ind


#def cluster_reads(splits, max_isize):    
#    trees = collections.defaultdict(lambda: ClusterTree(max_isize,1))
#    unmapped_reads = collections.defaultdict(lambda: [])
#    # bin the reads by reference and cluster by position
#    reads = []
#    read_ind = 0
#    for ind,split_reads in enumerate(splits):
#        ind_is_unmapped = True
#        for r in split_reads:
#            if r.is_unmapped:
#                continue
#            ind_is_unmapped = False
#            reads.append((ind,r))
#            trees[r.rname].insert(r.pos, r.aend, read_ind)
#            read_ind += 1
#        if ind_is_unmapped:
#            unmapped_reads[ind].extend(split_reads)
#    # find groups of reads on different references
#    for rname, cluster_tree in trees.iteritems():
#        ind_read_dict = collections.defaultdict(lambda: [])
#        for start, end, read_inds in cluster_tree.getregions():
#            for i in read_inds:
#                ind, r = reads[i]
#                ind_read_dict[ind].append(r)
#            ind_read_dict.update(unmapped_reads)
#            yield rname, start, end, ind_read_dict


def get_start_end_pos(mapped_inds, ind_read_dict):
    start = None
    end = None
    for ind in mapped_inds:        
        reads = ind_read_dict[ind]
        for r in reads:
            if start is None or r.pos < start:
                start = r.pos
            if end is None or r.aend > end:
                end = r.aend
    return start, end

def get_seq_and_qual(splits):
    seqs = []
    quals = []
    for reads in splits:
        # TODO: I am now manually switching genome hits to
        # 'unmapped' so that they are effectively ignored,
        # but they still need to be reverse complemented if 
        # seq is on the reverse strand, so using the 'unmapped'
        # flag alone does not suffice
        #if reads[0].is_unmapped:
        if not reads[0].is_reverse:
            seq = reads[0].seq
            qual = reads[0].qual
        else:
            # reverse complement sequences and quality scores 
            # back to original sequence read order
            seq = DNA_reverse_complement(reads[0].seq)
            qual = reads[0].qual[::-1]
        seqs.append(seq)
        quals.append(qual)        
    return ''.join(seqs), ''.join(quals)

#def get_fastq(partitions):
#    # all sequences are same so only use one partition
#    splits = partitions[0]
#    qname = splits[0][0].qname
#    seq, qual = get_seq_and_qual(splits)
#    return "@%s\n%s\n+%s\n%s" % (qname, seq, qname, qual)
#    #return ">%s\n%s" % (qname, seq)

def get_nonmapping_bedpe_string(pe_reads, read1_span, read2_span):
    qname = pe_reads[0][0][0][0].qname
    seq1, qual1 = get_seq_and_qual(pe_reads[0][0])
    seq2, qual2 = get_seq_and_qual(pe_reads[1][0])
    return '\t'.join('*', 0, 0, '*', 0, 0, qname, 1, 0, 0, seq1, seq2, 
                     qual1, qual2, int(read1_span), int(read2_span),
                     Chimera.NONMAPPING, 1) 


def make_chimera_mates(mate_splits, splits, gene_genome_map):
    '''
    returns a 3-tuple:
    * list of ChimeraMate objects,
    * whether the read has unmapped segments, and 
    * index of discordant splice
    '''
    refdict, mapped_inds, split_ind = find_first_split(splits)
    # read must have at least one mapped segment
    if len(mapped_inds) == 0:
        # this read has no mappings
        return [], True, split_ind
    # read has mappings, so extract information about the read    
    has_unmapped_splits = not set(xrange(len(mate_splits))).issubset(mapped_inds)
    span = mapped_inds[-1] < (len(mate_splits) - 1)
    seq, qual = get_seq_and_qual(mate_splits)    
    mates = []
    for rname, inddict in refdict.iteritems():
        g = gene_genome_map[rname]
        # TODO: there are cases where putative small rearrangements
        # create aberrant mappings of single reads to a gene, such
        # that it is plausible that the first segment aligns to a 
        # position larger than the last segment.  This is counterintuitive
        # and raises the need for indel and rearrangement calling prior to
        # fusion calling.  we mitigate this here by simply sorting the 
        # alignments in ascending order by position
        start, end = get_start_end_pos(mapped_inds, inddict)
        interval = Interval(start, end, chrom=g.tx_name, strand=g.strand)
        mates.append(ChimeraMate(interval, seq, qual, span))
    return mates, has_unmapped_splits, split_ind
    
def gen_chimera_candidates(splits5p, splits3p, gene_genome_map, read1_is_5prime):
    '''
    generator function yield (True, False, False, Chimera) objects if this is a 
    discordant pair.  otherwise, returns (False, spanning5p, spanning3p, None),
    where 'spanning5p' and 'spanning3p' are booleans indicating whether the
    5' and/or 3' splits should be realigned to detect spanning reads
    '''
    num_splits = len(splits5p) + len(splits3p)
    # find 5' genes and index of split where 3' genes begin
    mates5p, has_unmapped_splits, split_ind5p = make_chimera_mates(splits5p, 
                                                                   itertools.chain(splits5p, reversed(splits3p)),
                                                                   gene_genome_map)    
    concordant5p = False
    spanning5p = False
    if split_ind5p == num_splits:
        # no split found, so this is not a 
        # discordant read, but might be spanning
        concordant5p = True
        spanning5p = has_unmapped_splits
    # find 3' genes and index of split where 5' genes begin
    mates3p, has_unmapped_splits, split_ind3p = make_chimera_mates(splits3p,
                                                                   itertools.chain(splits3p, reversed(splits5p)), 
                                                                   gene_genome_map)
    concordant3p = False
    spanning3p = False
    if split_ind3p == num_splits:
        # no split found, so this is not a 
        # discordant read, but might be spanning
        concordant3p = True
        spanning3p = has_unmapped_splits
    if concordant5p or concordant3p:
        yield False, spanning5p, spanning3p, None
        return    
    # make chimera candidates
    qname = splits5p[0][0].qname    
    # determine the type of discordant read and also whether the reads could
    # be spanning the junction
    discordant_type = None
    if (split_ind5p < len(splits5p) and
        split_ind3p < len(splits3p)):
        discordant_type = Chimera.DISCORDANT_OVERLAPPING
    elif split_ind5p < len(splits5p):
        discordant_type = Chimera.DISCORDANT_5P
    elif split_ind3p < len(splits3p):
        discordant_type = Chimera.DISCORDANT_3P
    else:
        discordant_type = Chimera.DISCORDANT_INNER
    # produce all 5'/3' combinations
    for mate5p in mates5p:
        for mate3p in mates3p:
            yield True, False, False, Chimera(qname, discordant_type, mate5p, mate3p, read1_is_5prime)


def find_discordant_pairs(hits, gene_genome_map):
    '''
    function cannot be called with unmapped reads; at least one mapping
    must exist for each read in the pair
    '''
    # join hits
    chimeras = []
    for mate1, mate2 in ((0, 1), (1, 0)):
        mates5p = hits[0][mate1]
        mates3p = hits[1][mate2]
        for splits5p in mates5p:
            for splits3p in mates3p:                
                # combine 5'/3' partners
                for res in gen_chimera_candidates(splits5p, splits3p, 
                                                  gene_genome_map,
                                                  (mate1 == 0)):
                    is_discordant, span5p, span3p, chimera = res
                    if not is_discordant:
                        read1_span = (((mate1 == 0) and span5p) or
                                      ((mate1 == 1) and span3p))
                        read2_span = (((mate2 == 0) and span5p) or
                                      ((mate2 == 1) and span3p))                        
                        return [], read1_span, read2_span
                    chimeras.append(chimera)
    return chimeras, False, False

def discordant_reads_to_chimeras(input_bam_file, output_bedpe_file, gene_file,
                                 max_isize, library_type, contam_refs=None,
                                 spanning_reads_file=None):
    if contam_refs is None:
        contam_refs = []
    logging.info("Finding discordant reads")
    logging.debug("Input file: %s" % (input_bam_file))
    logging.debug("Output file: %s" % (output_bedpe_file))
    logging.debug("Spanning reads IDs file: %s" % (spanning_reads_file))
    logging.debug("Library type: %s" % (str(library_type)))
    logging.debug("Contaminant references: %s" % (contam_refs))
    # check the library type string is a tuple
    assert library_type == (0, 1)
    # build a map of gene name to genome coords
    logging.info("Reading gene index")
    bamfh = pysam.Samfile(input_bam_file, "rb")    
    gene_genome_map, gene_trees = build_gene_maps(bamfh, gene_file)
    # search for discordant reads that represent chimeras
    outfh = open(output_bedpe_file, "w")
    gene_tid_list = get_gene_tids(bamfh)
    # open unmapped fasta file if exists
    if spanning_reads_file is not None:
        spanningfh = open(spanning_reads_file, "w")
    else:
        spanningfh = None    
    # setup debugging logging messages
    debug_count = 0
    debug_every = 1e5
    debug_next = debug_every
    num_fragments = 0
    num_unmapped = 0
    num_unmapped_rescued = 0
    num_discordant = 0
    num_discordant_spanning = 0
    num_concordant = 0
    num_concordant_unmapped = 0
    logging.info("Parsing BAM file")    
    for pe_reads in parse_unpaired_reads(bamfh, gene_tid_list):
        # count fragments
        num_fragments += 1
        # progress log
        debug_count += 1
        if debug_count == debug_next:
            debug_next += debug_every
            logging.info("Total read pairs: %d" % (num_fragments))
            logging.info("Read pairs with at least one unmapped mate: %d (%d saved)" % (num_unmapped, num_unmapped_rescued))
            logging.info("Concordant reads: %d" % (num_concordant))            
            logging.info("Concordant pairs with at least one unmapped mate: %d" % (num_concordant_unmapped))
            logging.info("Discordant reads: %d" % (num_discordant))            
            logging.info("Discordant junction spanning reads: %d" % (num_discordant_spanning))            
        # bin reads by strand and mate and keep track of whether we
        # find mapped hits for this fragment
        hits = (([], []), ([], []))
        mate_mapping_codes = (set(), set())
        for mate, mate_partitions in enumerate(pe_reads):
            for strand, mapping_codes, splits in process_partition(mate_partitions):
                #print 'MATE', mate, 'STRAND', strand, 'CODES', mapping_codes
                #for i,s in enumerate(splits):
                #    for r in s:
                #        print 'INDEX', i, r
                mate_mapping_codes[mate].update(mapping_codes)
                if MAP in mapping_codes:
                    hits[strand][mate].append(splits)
        # if either mate is unmapped, stop here and do not proceed to
        # look for discordant pairs
        if any(MAP not in mapping_codes for mapping_codes in mate_mapping_codes):
            unmapped_mates = [False, False]            
            for mate, codes in enumerate(mate_mapping_codes):
                if NM in codes:
                    num_unmapped_rescued += 1
                    unmapped_mates[mate] = True
            # output unmapped reads for further testing (viruses, etc)
            # if the mate has any non-mapping reads (not genome/multimaps)
            print >>outfh, get_nonmapping_bedpe_string(pe_reads, 
                                                       unmapped_mates[0], 
                                                       unmapped_mates[1])            
            num_unmapped += 1
            continue
        # find discordant pairs in the mapped reads
        chimeras,any_read1_span,any_read2_span = find_discordant_pairs(hits, gene_genome_map)
        for chimera in chimeras:
            # output chimera
            print >>outfh, chimera.to_bedpe()
            # convert from 5'/3' back to read1/read2 for output
            if chimera.mate5p.is_spanning:
                any_read1_span = chimera.read1_is_5prime
                any_read2_span = not chimera.read1_is_5prime
            if chimera.mate3p.is_spanning:
                any_read1_span = not chimera.read1_is_5prime
                any_read2_span = chimera.read1_is_5prime
        # output putative spanning chimeras for further testing
        if (any_read1_span or any_read2_span):
            if len(chimeras) == 0:
                num_concordant_unmapped += 1
            else:             
                num_discordant_spanning += 1
            print >>outfh, get_nonmapping_bedpe_string(pe_reads, 
                                                       any_read1_span, 
                                                       any_read2_span) 
        if len(chimeras) > 0:
            num_discordant += 1
        else:
            num_concordant += 1
    # close output files
    outfh.close()
    if spanningfh is not None:
        spanningfh.close()
    # final progress
    logging.info("Total read pairs: %d" % (num_fragments))
    logging.info("Read pairs with at least one unmapped mate: %d (%d saved)" % (num_unmapped, num_unmapped_rescued))
    logging.info("Concordant reads: %d" % (num_concordant))            
    logging.info("Concordant pairs with at least one unmapped mate: %d" % (num_concordant_unmapped))
    logging.info("Discordant reads: %d" % (num_discordant))            
    logging.info("Discordant junction spanning reads: %d" % (num_discordant_spanning))  
    

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
    parser.add_option("--unmapped", dest="unmapped_reads_file", default=None)    
    options, args = parser.parse_args()
    input_bam_file = args[0]
    output_bedpe_file = args[1]
    gene_file = os.path.join(options.index_dir, config.GENE_FEATURE_FILE)
    library_type = parse_library_type(options.library_type)    
    discordant_reads_to_chimeras(input_bam_file, output_bedpe_file, gene_file, 
                                 options.max_fragment_length, library_type,
                                 options.contam_refs, options.unmapped_reads_file)

if __name__ == '__main__':
    main()



#
#def select_best_segments(seg_list):
#    # sort by number of unmapped segments
#    sorted_segs = []
#    for segs in seg_list:
#        num_unmapped_segs = sum(1 if r.is_unmapped else 0 for r in segs)
#        sorted_segs.append((num_unmapped_segs, segs))
#    sorted_segs = sorted(sorted_segs, key=operator.itemgetter(0))
#    best_num_unmapped_segs = sorted_segs[0][0]
#    best_segs = []
#    for num_unmapped_segs, segs in sorted_segs:
#        if num_unmapped_segs > best_num_unmapped_segs:
#            break
#        best_segs.append(segs)
#    return best_segs
