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

from bx.intersection import Interval, IntervalTree
import config
from base import parse_library_type
from feature import GeneFeature

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

def reorder_reads_5prime_to_3prime(partitions, tid_list):
    for partition in partitions:
        splits5p= collections.deque()
        splits3p = collections.deque()
        for split_reads in partition:
            reads5p = []
            reads3p = []
            # TODO: we select reads in the best 'strata', that is, the
            # set of reads with fewest mismatches to the reference.  is
            # this the best strategy, or should other metrics be employed
            # to choose from among multimapping reads?
            best_reads = select_best_mismatches(split_reads)
            for r in best_reads:
                # TODO: for now, we ignore genomic reads because they
                # require further processing
                if tid_list[r.rname] is None:
                    continue            
                # determine sense/antisense by assuming that
                # 5' reads are sense and 3' reads are antisense
                if r.is_unmapped:
                    reads3p.append(r)
                    reads5p.append(r)
                elif r.is_reverse:
                    reads3p.append(r)
                else:
                    reads5p.append(r)
            splits5p.append(reads5p)
            splits3p.appendleft(reads3p)        
        if all(len(reads) > 0 for reads in splits5p):
            yield 0, splits5p
        if all(len(reads) > 0 for reads in splits3p):
            yield 1, splits3p

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

ChimeraMate = collections.namedtuple('ChimeraMate',
                                     ['interval', 'seq', 'is_spanning'])

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
        mate5p = ChimeraMate(interval5p, fields[10], bool(int(fields[12])))
        interval3p = Interval(int(fields[4]),
                              int(fields[5]),
                              chrom=fields[3],
                              strand=fields[9])
        mate3p = ChimeraMate(interval3p, fields[11], bool(int(fields[13])))
        discordant_type = Chimera._discordant_types.index(fields[14])        
        read1_is_5prime = bool(int(fields[15]))        
        return Chimera(fields[6], discordant_type, mate5p, mate3p, read1_is_5prime)

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
    
def gen_chimera_candidates(splits5p, splits3p, gene_genome_map, read1_is_5prime):
    '''
    generator function yield (True, Chimera) objects if this is a 
    discordant pair and (False, None) otherwise
    '''
    splits = list(itertools.chain(splits5p, splits3p))
    # find 5' genes and index of split where 3' genes begin
    refdict5p, mapped_inds5p, split_ind5p = find_first_split(splits)
    if split_ind5p == len(splits):
        # no split found, so this is not a 
        # discordant read, but might be spanning
        yield False, None
    # lookup gene information for 5' partner
    span5p = mapped_inds5p[-1] < (len(splits5p) - 1)
    seq5p = ''.join(reads[0].seq for reads in splits5p)
    mates5p = []
    for rname, inddict in refdict5p.iteritems():        
        g = gene_genome_map[rname]
        # TODO: there are cases where putative small rearrangements
        # create aberrant mappings of single reads to a gene, such
        # that it is plausible that the first segment aligns to a 
        # position larger than the last segment.  This is counterintuitive
        # and raises the need for indel and rearrangement calling prior to
        # fusion calling.  we mitigate this here by simply sorting the 
        # alignments in ascending order by position
        start, end = get_start_end_pos(mapped_inds5p, inddict)
        interval = Interval(start, end, chrom=g.tx_name, strand=g.strand)
        mates5p.append(ChimeraMate(interval, seq5p, span5p))
    # find 3' genes and index of split where 5' genes begin
    refdict3p, mapped_inds3p, split_ind3p = find_first_split(reversed(splits))
    # lookup gene information for 3' partner
    span3p = mapped_inds3p[-1] < (len(splits3p) - 1)
    seq3p = ''.join(reads[0].seq for reads in splits3p)
    mates3p = []
    for rname, inddict in refdict3p.iteritems():
        g = gene_genome_map[rname]
        # 3' is antisense, so 'last' index should be closest to start
        # and 'first' index should be closest to end of transcript
        # TODO: we use the sorting scheme here rather than intuition because
        # some corner cases violate the principles (see above)
        start, end = get_start_end_pos(mapped_inds3p, inddict)
        interval = Interval(start, end, chrom=g.tx_name, strand=g.strand)
        mates3p.append(ChimeraMate(interval, seq3p, span3p))
    #
    # make chimera candidates
    #
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
            yield True, Chimera(qname, discordant_type, mate5p, mate3p, read1_is_5prime)


def check_read_unmapped(split_partitions):
    # if only one partition with one split
    # and one mapping read that is unmapped,
    # then the entire read is unmapped
    if (len(split_partitions) == 1 and
        len(split_partitions[0]) == 1 and
        len(split_partitions[0][0]) == 1 and
        split_partitions[0][0][0].is_unmapped):
        r = split_partitions[0][0][0]
        multimap = r.opt('XM') > 0
        return True, multimap
    else:
        return False, False

def get_fasta(partitions):
    # TODO: reverse complement reverse sequences to 
    # ensure quality scores are in original sequenced
    # read order?
    #from seq import DNA_reverse_complement
    # all sequences are same so only use one partition
    splits = partitions[0]
    qname = splits[0][0].qname
    seq = ''.join(reads[0].seq for reads in splits)
    return ">%s\n%s" % (qname, seq)
    #qual = ''.join(reads[0].qual for reads in splits)
    #return "@%s\n%s\n+%s\n%s" % (qname, seq, qname, qual)

def find_discordant_pairs(pe_reads, tid_list, genome_tid_set, 
                          gene_genome_map, gene_trees, max_isize, 
                          contam_tids):
    '''
    function cannot be called with unmapped reads; at least one mapping
    must exist for each read in the pair
    '''
    # reorder reads in 5' -> 3' direction and bin by strand and mate
    hits = (([], []), ([], []))
    for mate, mate_partitions in enumerate(pe_reads):
        for strand,splits in reorder_reads_5prime_to_3prime(mate_partitions, tid_list):
            hits[strand][mate].append(splits)
    # join hits
    chimeras = []
    for mate1, mate2 in ((0, 1), (1, 0)):
        for splits5p in hits[0][mate1]:
            for splits3p in hits[1][mate2]:
                # combined 5'/3' partners
                for is_discordant,chimera in gen_chimera_candidates(splits5p, 
                                                                    splits3p, 
                                                                    gene_genome_map,
                                                                    (mate1 == 0)):
                    if not is_discordant:
                        return []
                    chimeras.append(chimera)
    return chimeras

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
                                 unmapped_fasta_file=None):
    if contam_refs is None:
        contam_refs = []
    logging.info("Finding discordant reads")
    logging.debug("Input file: %s" % (input_bam_file))
    logging.debug("Output file: %s" % (output_bedpe_file))
    logging.debug("Unmapped reads FASTA file: %s" % (unmapped_fasta_file))
    logging.debug("Library type: %s" % (library_type))
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
    if unmapped_fasta_file is not None:
        fastafh = open(unmapped_fasta_file, "w")
    else:
        fastafh = None    
    # setup debugging logging messages
    debug_count = 0
    debug_every = 1e5
    debug_next = debug_every
    num_fragments = 0
    num_unmapped = 0
    num_discordant = 0
    num_discordant_spanning = 0
    num_concordant = 0
    for pe_reads in parse_unpaired_reads(bamfh):
        # check that reads map
        any_is_unmapped = False
        for mate, mate_partitions in enumerate(pe_reads):        
            is_unmapped, is_multimap = check_read_unmapped(mate_partitions)
            any_is_unmapped = any_is_unmapped or is_unmapped
            # output unmapped reads for further testing (viruses, etc)
            if (fastafh is not None) and is_unmapped and not is_multimap:
                print >>fastafh, get_fasta(mate_partitions) 
        if any_is_unmapped:
            num_unmapped += 1
        else:
            chimeras = find_discordant_pairs(pe_reads, gene_tid_list, genome_tid_set, gene_genome_map, 
                                             gene_trees, max_isize, contam_tids)
            any_read1_span = False
            any_read2_span = False            
            for chimera in chimeras:
                print >>outfh, chimera.to_bedpe()
                if chimera.mate5p.is_spanning:
                    if chimera.read1_is_5prime:
                        any_read1_span = True
                    else:
                        any_read2_span = True
                if chimera.mate3p.is_spanning:
                    if chimera.read1_is_5prime:
                        any_read1_span = True
                    else:
                        any_read2_span = True
            # output putative spanning chimeras for further testing
            if (fastafh is not None) and (any_read1_span or any_read2_span):                
                num_discordant_spanning += 1
                if any_read1_span:
                    print >>fastafh, ">%s\n%s" % (chimeras[0].qname, chimeras[0].mate5p.seq)
                if any_read2_span:
                    print >>fastafh, ">%s\n%s" % (chimeras[0].qname, chimeras[0].mate3p.seq)
            if len(chimeras) > 0:
                num_discordant += 1
            else:
                num_concordant += 1
        num_fragments += 1
        # progress log
        debug_count += 1
        if debug_count == debug_next:
            debug_next += debug_every
            logging.info("Total read pairs: %d" % (num_fragments))
            logging.info("Read pairs with at least one unmapped mate: %d" % (num_unmapped))
            logging.info("Discordant reads: %d" % (num_discordant))            
            logging.info("Discordant junction spanning reads: %d" % (num_discordant_spanning))            
            logging.info("Concordant reads: %d" % (num_concordant))            
    # close output files
    outfh.close()
    if fastafh is not None:
        fastafh.close()
    # final progress
    logging.info("Total read pairs: %d" % (num_fragments))
    logging.info("Read pairs with at least one unmapped mate: %d" % (num_unmapped))
    logging.info("Discordant reads: %d" % (num_discordant))           
    logging.info("Discordant junction spanning reads: %d" % (num_discordant_spanning))            
    logging.info("Concordant reads: %d" % (num_concordant))            


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
    parser.add_option("--unmapped", dest="unmapped_fasta_file", default=None)    
    options, args = parser.parse_args()
    input_bam_file = args[0]
    output_bedpe_file = args[1]
    gene_file = os.path.join(options.index_dir, config.GENE_FEATURE_FILE)
    library_type = parse_library_type(options.library_type)    
    discordant_reads_to_chimeras(input_bam_file, output_bedpe_file, gene_file, 
                                 options.max_fragment_length, options.library_type,
                                 options.contam_refs, options.unmapped_fasta_file)

if __name__ == '__main__':
    main()

#def find_discordant_pairs(pe_reads, tid_list, genome_tid_set, 
#                          gene_genome_map, gene_trees, max_isize, 
#                          contam_tids):
#    filtered_pe_reads = filter_genomic_reads(pe_reads, tid_list)
#    # find out whether the reads by themselves are concordant first
#    r1_strand_rnames, r1_unmapped, r1_concordant, r1_discordant = get_discordant_partitions(filtered_pe_reads[0], tid_list)
#    r2_strand_rnames, r2_unmapped, r2_concordant, r2_discordant = get_discordant_partitions(filtered_pe_reads[1], tid_list)
#    #r1_is_concordant = len(r1_concordant) > 0
#    #r2_is_concordant = len(r2_concordant) > 0
#    r1_is_discordant = len(r1_concordant) == 0 and len(r1_discordant) > 0
#    r2_is_discordant = len(r2_concordant) == 0 and len(r2_discordant) > 0
#    #print 'r1', r1_rnames, r1_nm, r1_concordant, r1_discordant
#    #print 'r2', r2_rnames, r2_nm, r2_concordant, r2_discordant
#    if len(r1_unmapped) > 0 or len(r2_unmapped) > 0:        
#        # if either read is nonmapping cannot predict discordant
#        print 'ONE OR BOTH UNMAPPED'
#        print 'ONE OR BOTH UNMAPPED'
#        print 'ONE OR BOTH UNMAPPED'
#        return    
#    elif (not r1_is_discordant) and (not r2_is_discordant):
#        # both reads are concordant, so check their reference names
#        # for matches on opposite strands
#        is_concordant = False
#        for strands in ((0,1), (1,0)):
#            shared_rnames = r1_strand_rnames[strands[0]].intersection(r2_strand_rnames[strands[1]])
#            if len(shared_rnames) > 0:
#                is_concordant = True
#                break
#        if is_concordant:
#            # reads are concordant
#            print 'READS CONCORDANT'
#            print 'READS CONCORDANT'
#            print 'READS CONCORDANT'
#            return    
#    print 'READS DISCORDANT'
#    print 'READS DISCORDANT'
#    print 'READS DISCORDANT'
#    # this is a discordant pair
#
#    # reorder reads in 5' -> 3' direction and bin by strand and mate
#    hits = (([], []), ([], []))
#    for mate, mate_partitions in enumerate(filtered_pe_reads):
#        for strand,splits in reorder_reads_5prime_to_3prime(mate_partitions):
#            hits[strand][mate].append(splits)
#
#    # join hits
#    for mate1, mate2 in ((0, 1), (1, 0)):
#        for splits5p in hits[0][mate1]:
#            for splits3p in hits[1][mate2]:
#                # combined 5'/3' partners
#                make_chimera_candidate(splits5p, splits3p)

#    # select hits with fewest mismatches, keeping categories separate for now
#    concordant_partially_mapped = select_best_mismatches(concordant_partially_mapped)
#    discordant_unmapped = select_best_mismatches(discordant_unmapped)
#    discordant_mapped = select_best_mismatches(discordant_mapped)
#    
#    #return discordant_mapped, discordant_unmapped
#
#    for segs in discordant_mapped:
#        print '|||'.join(['%s:%d' % ('None' if r.rname == -1 else tid_list[r.rname], r.pos) for r in segs])
#
#    #print 'mate', mate1, '|||'.join(['%s:%d' % (tid_list[r.rname], r.pos) for r in fiveprime_segs])
#    #print 'mate', mate2, '|||'.join(['%s:%d' % (tid_list[r.rname], r.pos) for r in threeprime_segs])
#    #print fiveprime_segs
#    #print threeprime_segs
#
#    for mate_hits in pe_reads:
#        for split_partitions in mate_hits:
#            for split_reads in split_partitions:
#                for r in split_reads:
#                    # ignore unmapped reads
#                    if r.is_unmapped:
#                        continue
#                    # ignore genomic reads
#                    if tid_list[r.rname] is None:
#                        continue
#                    output_bamfh.write(r)

#def get_strand_and_ref_names(partition_reads, tid_list):    
#    strand_rnames = ([],[])
#    for split_reads in partition_reads:        
#        split_strand_rnames = (set(), set())
#        for r in split_reads:
#            if r.is_unmapped:
#                continue
#            split_strand_rnames[int(r.is_reverse)].add(r.rname)
#        strand_rnames[0].append(split_strand_rnames[0])
#        strand_rnames[1].append(split_strand_rnames[1])
#    return strand_rnames
#
#def check_strand_concordant(strand_rnames):
#    shared_rnames = set()
#    split_ind = 0
#    for split_rnames in strand_rnames:
#        split_ind += 1
#        if len(split_rnames) > 0:
#            shared_rnames.update(split_rnames)
#            break
#    for split_rnames in strand_rnames[split_ind:]:
#        if len(split_rnames) > 0:
#            shared_rnames.intersection_update(split_rnames)
#            if len(shared_rnames) == 0:
#                return False, split_ind
#        split_ind += 1
#    return True, len(strand_rnames)
#
#def get_discordant_partitions(read_hits, tid_list):
#    discordant_partitions = []
#    concordant_partitions = []
#    nonmapping_partitions = []
#    all_strand_rnames = (set(), set())
#    for partition_ind, partition_reads in enumerate(read_hits):        
#        partition_strand_rnames = get_strand_and_ref_names(partition_reads, tid_list)        
#        all_strand_rnames[0].update(*partition_strand_rnames[0])
#        all_strand_rnames[1].update(*partition_strand_rnames[1])        
#        # if there are no reference names then this read is unmapped
#        if all(len(x) == 0 for x in all_strand_rnames):
#            nonmapping_partitions.append((partition_ind, 0, 0))
#            continue        
#        for strand, strand_rnames in enumerate(partition_strand_rnames):
#            is_concordant, split_ind = check_strand_concordant(strand_rnames)
#            if is_concordant:                
#                concordant_partitions.append((partition_ind, strand, split_ind))
#            else:
#                discordant_partitions.append((partition_ind, strand, split_ind))
#    return all_strand_rnames, nonmapping_partitions, concordant_partitions, discordant_partitions
#
#
#def filter_genomic_reads(pe_reads, tid_list):    
#    new_pe_reads = ([],[])
#    for mate,mate_partitions in enumerate(pe_reads):
#        for partition in mate_partitions:
#            new_partition = []
#            for split_reads in partition:
#                new_split_reads = [r for r in split_reads 
#                                   if tid_list[r.rname] is not None]
#                new_partition.append(new_split_reads)
#            new_pe_reads[mate].append(new_partition)
#    return new_pe_reads
