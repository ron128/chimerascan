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
    
    def update_missing(self, other_dict):
        for ind, reads in other_dict.iteritems():
            if ind not in self.split_read_dict:
                self.split_read_dict[ind] = reads
                 
    @staticmethod
    def get_nonmapping(split_read_dict):
        yield ReadCluster(-1, 0, 0, 0, split_read_dict)

class RefMap(object):
    def __init__(self, rname, max_dist):
        self.rname = rname
        self.cluster_trees = (ClusterTree(max_dist,1),
                              ClusterTree(max_dist,1))                             
        self.reads = []

    def add_read(self, r):
        strand = int(r.is_reverse)
        self.cluster_trees[strand].insert(r.pos, r.aend, len(self.reads))
        self.reads.append(r)

    def get_read_clusters(self):
        for strand, cluster_tree in enumerate(self.cluster_trees):
            for start, end, read_inds in cluster_tree.getregions():
                split_read_dict = collections.defaultdict(lambda: [])
                # organize reads by split index
                for i in read_inds:
                    r = self.reads[i]
                    split_read_dict[r.opt(RTAG_SPLIT_IND)].append(r)
                if len(split_read_dict) > 0:
                    yield ReadCluster(self.rname, start, end, strand, split_read_dict)


class DiscordantPair(object):
    def __init__(self):
        pass

        
def find_split_read_clusters(partition_splits, max_indel_size, max_multihits):
    refmaps = {}
    split_mapping_codes = []
    unmapped_read_dict = collections.defaultdict(lambda: [])
    for split_ind, split_reads in enumerate(partition_splits):
        mapping_codes = set()                         
        for r in split_reads:
            # keep track of mapping results for reads in this split
            #print 'IND', split_ind, 'mapping code', _mapping_code_strings[get_mapping_code(r, max_multihits)]            
            mapping_codes.add(get_mapping_code(r, max_multihits))
            if r.is_unmapped:
                unmapped_read_dict[split_ind].append(r)
                continue
            # cluster reads by reference names
            if r.rname not in refmaps:
                refmaps[r.rname] = RefMap(r.rname, max_indel_size)
            refmaps[r.rname].add_read(r)
        split_mapping_codes.append(mapping_codes)
    # search read clusters to find concordant clusters
    rclusters = []
    ind_clust_map = collections.defaultdict(lambda: set())
    for rname, refmap in refmaps.iteritems():
        for rclust in refmap.get_read_clusters():
            # add to master list of read clusters
            clust_id = len(rclusters)
            rclusters.append(rclust)
            # build an index from split index to list of read clusters
            # that contain that index 
            for ind in rclust.split_read_dict:
                ind_clust_map[ind].add(clust_id)
            # fill in missing splits in cluster with unmapped reads
            # AFTER adding indexes to map so that they do not affect
            # detection of discordant reads
            rclust.update_missing(unmapped_read_dict)
    # now walk through split indexes in order and find split points
    start_ind = 0
    last_mapped_ind = 0
    current_clusters = set()
    concordant_clust_ids = []
    for split_ind in xrange(len(partition_splits)):
        # if there are no clusters then this index
        # is unmapped and not useful
        if split_ind not in ind_clust_map:
            continue  
        # get read clusters at this index
        clust_ids = set(ind_clust_map[split_ind])
        if len(current_clusters) == 0:
            # initialize the current clusters to the 
            # first mapped index
            current_clusters.update(clust_ids)
        elif current_clusters.isdisjoint(clust_ids):
            # save the set of concordant clusters and the split index
            concordant_clust_ids.append((start_ind, split_ind, current_clusters))
            # initialize new clusters
            current_clusters = clust_ids            
            start_ind = last_mapped_ind + 1
        else:
            # intersect index clusters together to 
            # reduce number of current clusters
            current_clusters.intersection_update(clust_ids)
        # keep track of the last mapped index so that any
        # unmapped indexes in between can be attributed to 
        # both the 5' and 3' split genes
        last_mapped_ind = split_ind
    # cleanup the last set of clusters
    if len(current_clusters) > 0:
        concordant_clust_ids.append((start_ind, split_ind+1, current_clusters))
    # convert from clusters to read segments
    concordant_clusters = []
    for start_ind, end_ind, clust_ids in concordant_clust_ids:
        concordant_clusters.append((start_ind, end_ind, tuple(rclusters[id] for id in clust_ids)))
#        print 'START', start_ind, 'END', end_ind, 'IDS', concordant_clust_ids
#        print 'LEN', len(rclusters)
#        for id in clust_ids:
#            c = rclusters[id]
#            print 'CLUSTER'
#            for i in xrange(start_ind, end_ind):
#                print c.split_read_dict[i]
    # fill in unmapped reads within cluster
    #rclust.update_missing(unmapped_read_dict)        
    #print 'UNMAPPED', dict(unmapped_read_dict)
    #print 'IND', split_ind, split_mapping_codes[split_ind], ind_clust_map[split_ind]
    return concordant_clusters


def bin_alignments_by_reference(pe_reads, max_indel_size, max_isize, max_multihits):
    for mate, partitions in enumerate(pe_reads):
        print 'MATE', mate
        for partition_ind, partition_splits in enumerate(partitions):
            print 'PARTITION', partition_ind
            print 'NUM SPLITS', len(partition_splits)            
            find_split_read_clusters(partition_splits, max_isize, max_multihits)


def find_discordant_reads(bamfh, max_indel_size, max_isize, max_multihits):
    for pe_reads in parse_reads(bamfh):
        bin_alignments_by_reference(pe_reads, max_indel_size, max_isize, max_multihits)

def interval_overlap(chrom1, start1, end1, chrom2, start2, end2):
    if chrom1 != chrom2:
        return False
    return (start1 < end2) and (start2 < end1)    


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
