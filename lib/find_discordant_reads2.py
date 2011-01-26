'''
Created on Jan 25, 2011

@author: mkiyer
'''
import collections
import logging
import operator
import os

# local libs
import pysam
from bx.intersection import Interval, IntervalTree
from bx.cluster import ClusterTree

# local imports
import config
from base import parse_library_type
from seq import DNA_reverse_complement
from gene_to_genome import build_gene_maps, get_gene_tids


def parse_reads(bamfh):
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


# Mapping codes
NONMAPPING = 0
MULTIMAPPING = 1
MAPPING = 2
_mapping_code_strings = ["NONMAPPING", "MULTIMAPPING", "MAPPING"]

# Read tags
RTAG_NUM_PARTITIONS = "NH"
RTAG_PARTITION_IND = "XH"
RTAG_NUM_SPLITS = "XN"
RTAG_SPLIT_IND = "XI"
RTAG_NUM_MAPPINGS = "IH"
RTAG_MAPPING_IND = "HI"

def get_mapping_code(read, multihit_limit):
    if read.is_unmapped:
        if read.opt('XM') >= multihit_limit:
            return MULTIMAPPING
        else:
            return NONMAPPING
    return MAPPING

class ReadMetadata(object):
    __slots__ = ('read', 'mapping_code')
    
    def __init__(self, r, multihit_limit=1):
        self.read = r
        self.mapping_code = get_mapping_code(r, multihit_limit)

class ReadCluster(object):
    def __init__(self, rname, start, end, strand, split_read_dict):
        self.rname = rname
        self.start = start
        self.end = end
        self.strand = strand
        self.split_read_dict = split_read_dict
    
    @staticmethod
    def get_nonmapping(split_read_dict):
        yield ReadCluster(-1, 0, 0, 0, split_read_dict)


def interval_overlap(chrom1, start1, end1, chrom2, start2, end2):
    if chrom1 != chrom2:
        return False
    return (start1 < end2) and (start2 < end1)    

class RefMap(object):
    def __init__(self, rname, max_dist):
        self.rname = rname
        self.cluster_tree = ClusterTree(max_dist,1)
        self.reads = []

    def add_read(self, r):
        self.cluster_tree.insert(r.pos, r.aend, len(self.reads))
        self.reads.append(r)

    def get_read_clusters(self):
        for start, end, read_inds in self.cluster_tree.getregions():
            strand_split_read_dicts = (collections.defaultdict(lambda: []),
                                       collections.defaultdict(lambda: []))
            # organize reads by split index
            for i in read_inds:
                r = self.reads[i]
                strand = int(r.is_reverse)
                strand_split_read_dicts[strand][r.opt(RTAG_SPLIT_IND)].append(r)
            for strand, split_read_dict in strand_split_read_dicts:
                if len(split_read_dict) > 0:
                    yield ReadCluster(self.rname, start, end, strand, split_read_dict)


class DiscordantPair(object):
    def __init__(self):
        pass

        
def bin_partition_by_reference(partition_splits, max_isize):
    refmaps = {}
    split_mapping_codes = []
    unmapped_read_dict = collections.defaultdict(lambda: [])
    for split_ind, split_reads in enumerate(partition_splits):
        mapping_codes = set()                         
        for r in split_reads:
            # keep track of mapping results for reads in this split
            mapping_codes.add(get_mapping_code(r))
            if r.is_unmapped:
                unmapped_read_dict[split_ind].append(r)
                continue
            # cluster reads by reference names
            if r.rname not in refmaps:
                refmaps[r.rname] = RefMap(r.rname, max_isize)
            refmaps[r.rname].add_read(r)
        split_mapping_codes.append(mapping_codes)
    for rname, refmap in refmaps.iteritems():
        rclusters = list(RefMap.get_read_clusters())
    

def bin_alignments_by_reference(pe_reads, max_isize):
    for mate, partitions in enumerate(pe_reads):
        for partition_ind, partition_splits in enumerate(partitions):
            bin_partition_by_reference(partition_splits, max_isize)


def find_discordant_reads(bamfh, max_indel_size, max_isize, max_multihits):
    for pe_reads in parse_reads(bamfh):
        print pe_reads

def main():
    from optparse import OptionParser
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = OptionParser("usage: %prog [options] <bam> <out.bedpe>")
    parser.add_option('--max-fragment-length', dest="max_fragment_length", 
                      type="int", default=1000)
    parser.add_option('--max-indel-size', dest="max_indel_size", 
                      type="int", default=100)
    parser.add_option('--library-type', dest="library_type", default="fr")
    parser.add_option("--index", dest="index_dir",
                      help="Path to chimerascan index directory")
    parser.add_option('--multihits', type="int", default=1)
    parser.add_option("--contam-refs", dest="contam_refs", default=None)
    options, args = parser.parse_args()
    input_bam_file = args[0]
    output_bedpe_file = args[1]
    gene_file = os.path.join(options.index_dir, config.GENE_FEATURE_FILE)
    library_type = parse_library_type(options.library_type)    

    bamfh = pysam.Samfile(input_bam_file, "rb")
    find_discordant_reads(bamfh, 
                          max_indel_size=options.max_indel_size, 
                          max_isize=options.max_fragment_length,
                          max_multihits=options.multihits)                          
    bamfh.close()

if __name__ == '__main__':
    main()



