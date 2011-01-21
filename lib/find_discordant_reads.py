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

import pysam

# local imports
from bx.intersection import Interval, IntervalTree
from bx.cluster import ClusterTree

import config
from base import parse_library_type
from feature import GeneFeature
from seq import DNA_reverse_complement

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
    _discordant_types = ["DISCORDANT_INNER",
                         "DISCORDANT_5PRIME",
                         "DISCORDANT_3PRIME",
                         "DISCORDANT_OVERLAPPING"]
    
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

def parse_unpaired_reads(bamfh):
    pe_reads = ([], [])
    # reads must be binned by qname, mate, hit, and segment
    # so initialize to mate 0, hit 0, segment 0
    num_reads = 0
    prev_qname = None
    for read in bamfh:
        # ignore paired reads
        if read.is_proper_pair:
            continue
        # get read attributes
        qname = read.qname
        mate = 0 if read.is_read1 else 1
        # get hit/segment/mapping tags
        num_split_partitions = read.opt('NH')
        partition_ind = read.opt('XH')
        num_splits = read.opt('XN')
        split_ind = read.opt('XI')
        num_mappings = read.opt('IH')
        mapping_ind = read.opt('HI')        
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
            pe_reads[mate].extend([list() for x in xrange(num_split_partitions)])
        mate_reads = pe_reads[mate]
        # initialize hit segments
        if len(mate_reads[partition_ind]) == 0:
            mate_reads[partition_ind].extend([list() for x in xrange(num_splits)])
        split_reads = mate_reads[partition_ind][split_ind]
        # initialize segment mappings
        if len(split_reads) == 0:
            split_reads.extend([None for x in xrange(num_mappings)])
        # add segment to hit/mate/read
        split_reads[mapping_ind] = read
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

def process_partition(partitions, tid_list):
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
                # TODO: for now, we treat genomic reads as unmapped because they
                # require further processing.  by treating unmapped we allow them
                # to be consider as mis-mapped spanning reads in future steps,
                # and allow the other segments in the paired-end fragment to 
                # determine  
                if tid_list[r.rname] is None:
                    r.is_unmapped = True
                    reads5p.append(r)
                    reads3p.append(r)
                    codes5p.add(GENOME)
                    codes3p.add(GENOME)
                elif r.is_unmapped:
                    reads3p.append(r)
                    reads5p.append(r)
                    code = MULTIMAP if r.opt('XM') > 0 else NM
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
            print 'INDEX', ind, 'READ', r
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

def get_fastq(partitions):
    # all sequences are same so only use one partition
    splits = partitions[0]
    qname = splits[0][0].qname
    seq, qual = get_seq_and_qual(splits)
    return "@%s\n%s\n+%s\n%s" % (qname, seq, qname, qual)
    #return ">%s\n%s" % (qname, seq)

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
    mates = []
    span = mapped_inds[-1] < (len(mate_splits) - 1)
    seq, qual = get_seq_and_qual(mate_splits)     
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
    return mates, mapped_inds[0] > 0, split_ind
    
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
    print 'NUM SPLITS', num_splits
    print 'MATES5P', mates5p
    print 'UNMAPPED', has_unmapped_splits
    print '5P SPLIT IND', split_ind5p
    
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
    print 'MATES3P', mates3p
    print 'UNMAPPED', has_unmapped_splits
    print '3P SPLIT IND', split_ind3p

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

def discordant_reads_to_chimeras(input_bam_file, output_bedpe_file, gene_file,
                                 max_isize, library_type, contam_refs=None,
                                 unmapped_fastq_file=None):
    if contam_refs is None:
        contam_refs = []
    logging.info("Finding discordant reads")
    logging.debug("Input file: %s" % (input_bam_file))
    logging.debug("Output file: %s" % (output_bedpe_file))
    logging.debug("Unmapped reads FASTA file: %s" % (unmapped_fastq_file))
    logging.debug("Library type: %s" % (str(library_type)))
    logging.debug("Contaminant references: %s" % (contam_refs))
    # build a map of gene name to genome coords
    logging.info("Reading gene index")
    bamfh = pysam.Samfile(input_bam_file, "rb")    
    gene_genome_map, gene_trees = build_gene_maps(bamfh, gene_file)
    # check the library type string into a tuple
    assert library_type == (0, 1)
    #same_strand = (library_type[0] == library_type[1])    
    # get contaminant reference tids
    if contam_refs is None:
        contam_tids = set()
    else:
        contam_tids = set(get_tids(bamfh, contam_refs))
    # search for discordant reads that represent chimeras
    outfh = open(output_bedpe_file, "w")
    gene_tid_list = get_gene_tids(bamfh)
    genome_tid_set = get_genome_tids(bamfh)
    # open unmapped fasta file if exists
    if unmapped_fastq_file is not None:
        fastqfh = open(unmapped_fastq_file, "w")
    else:
        fastqfh = None    
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
    for pe_reads in parse_unpaired_reads(bamfh):
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
            for strand, mapping_codes, splits in process_partition(mate_partitions, gene_tid_list):
#                print 'MATE', mate, 'STRAND', strand, 'CODES', mapping_codes
#                for i,s in enumerate(splits):
#                    for r in s:
#                        print 'INDEX', i, r
                mate_mapping_codes[mate].update(mapping_codes)
                if MAP in mapping_codes:
                    hits[strand][mate].append(splits)
        # if either mate is unmapped, stop here and do not proceed to
        # look for discordant pairs
        if any(MAP not in mapping_codes for mapping_codes in mate_mapping_codes):
            for mate, codes in enumerate(mate_mapping_codes):
                # output unmapped reads for further testing (viruses, etc)
                # if the mate has any non-mapping reads (not genome/multimaps)
                if NM in codes:
                    num_unmapped_rescued += 1
                    if (fastqfh is not None):                    
                        print >>fastqfh, get_fastq(pe_reads[mate])
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
        if (fastqfh is not None) and (any_read1_span or any_read2_span):
            if len(chimeras) == 0:
                num_concordant_unmapped += 1
            else:             
                num_discordant_spanning += 1
            if any_read1_span:
                print >>fastqfh, get_fastq(pe_reads[0]) 
            if any_read2_span:
                print >>fastqfh, get_fastq(pe_reads[1])
        if len(chimeras) > 0:
            num_discordant += 1
        else:
            num_concordant += 1
    # close output files
    outfh.close()
    if fastqfh is not None:
        fastqfh.close()
    # final progress
    logging.info("Total read pairs: %d" % (num_fragments))
    logging.info("Read pairs with at least one unmapped mate: %d" % (num_unmapped))
    logging.info("Discordant reads: %d" % (num_discordant))           
    logging.info("Discordant junction spanning reads: %d" % (num_discordant_spanning))            
    logging.info("Concordant reads: %d" % (num_concordant))            
    logging.info("Concordant pairs with at least one unmapped mate: %d" % (num_concordant_unmapped))            


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
    parser.add_option("--unmapped", dest="unmapped_fastq_file", default=None)    
    options, args = parser.parse_args()
    input_bam_file = args[0]
    output_bedpe_file = args[1]
    gene_file = os.path.join(options.index_dir, config.GENE_FEATURE_FILE)
    library_type = parse_library_type(options.library_type)    
    discordant_reads_to_chimeras(input_bam_file, output_bedpe_file, gene_file, 
                                 options.max_fragment_length, library_type,
                                 options.contam_refs, options.unmapped_fastq_file)

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
