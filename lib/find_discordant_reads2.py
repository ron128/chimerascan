'''
Created on Jan 25, 2011

@author: mkiyer
'''
import collections
import itertools
import logging
import operator
import os

# local libs
import pysam
from bx.intersection import Interval, IntervalTree
from bx.cluster import ClusterTree

# local imports
import config
from base import parse_library_type, SamTags
from gene_to_genome import build_gene_maps

# Mapping codes
NONMAPPING = 0
MULTIMAPPING = 1
MAPPING = 2
_mapping_code_strings = ["NONMAPPING", "MULTIMAPPING", "MAPPING"]
    
def get_mapping_code(read, multihit_limit):
    if read.is_unmapped:
        if read.opt(SamTags.RTAG_BOWTIE_MULTIMAP) >= multihit_limit:
            return MULTIMAPPING
        else:
            return NONMAPPING
    return MAPPING

def parse_reads(bamfh): 
    # reads must be binned by qname, mate, hit, and segment
    # so initialize to mate 0, hit 0, segment 0
    pe_reads = ([], [])
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
        num_split_partitions = read.opt(SamTags.RTAG_NUM_PARTITIONS)
        partition_ind = read.opt(SamTags.RTAG_PARTITION_IND)
        num_splits = read.opt(SamTags.RTAG_NUM_SPLITS)
        split_ind = read.opt(SamTags.RTAG_SPLIT_IND)
        num_mappings = read.opt(SamTags.RTAG_NUM_MAPPINGS)
        mapping_ind = read.opt(SamTags.RTAG_MAPPING_IND)
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



class ReadMetadata(object):
    __slots__ = ('read', 'mapping_code')
    
    def __init__(self, r, multihit_limit=1):
        self.read = r
        self.mapping_code = get_mapping_code(r, multihit_limit)

class RefCluster(object):
    def __init__(self, rname, strand, max_dist):
        self.rname = rname
        self.strand = strand
        self.cluster_tree = ClusterTree(max_dist,1)
        self.vals = []

    def add(self, start, end, ind, val):
        self.cluster_tree.insert(start, end, len(self.vals))
        self.vals.append((val, ind))

    def get_clusters(self):
        for start, end, inds in self.cluster_tree.getregions():
            ind_val_dict = collections.defaultdict(lambda: [])
            # organize reads by split index
            for i in inds:
                val,ind = self.vals[i]
                ind_val_dict[ind].append(val)
            # TODO: remove assert
            assert len(ind_val_dict) > 0
            yield Interval(start, end, chrom=self.rname, strand=self.strand, value=ind_val_dict)

class ReadCluster(object):
    def __init__(self, interval, start_ind, end_ind, mapped_inds, 
                 unmapped_inds, unmapped_read_dict):
        # access interval object with start, end, strand, chrom info
        self.start = interval.start
        self.end = interval.end
        self.rname = interval.chrom
        self.strand = interval.strand
        self.mate = None
        # minimum multimaps for all reads in cluster
        self.multimaps = None
        self.split_dict = interval.value
        for ind,reads in interval.value.iteritems():
            min_multimaps = min(r.opt('NH') for r in reads)
            if (self.multimaps is None) or self.multimaps < min_multimaps:
                self.multimaps = min_multimaps
        # Indexes within split segments in the cluster
        self.start_ind = start_ind
        self.end_ind = end_ind
        self.mapped_inds = mapped_inds
        self.unmapped_inds = unmapped_inds
        # fill in missing splits in cluster with unmapped splits
        for ind in unmapped_inds:
            # TODO: remove assert
            assert ind not in self.split_dict
            self.split_dict[ind] = unmapped_read_dict[ind]

    def __repr__(self):
        return ("<%s(rname=%s, strand=%s, start=%d, end=%d, multimaps=%d, mate=%d, start_ind=%d, "
                "end_ind=%d, mapped_inds=%s, unmapped_inds=%s)>" % 
                (self.__class__.__name__, self.rname, self.strand, self.start, self.end,
                 self.multimaps, self.mate, self.start_ind, self.end_ind, self.mapped_inds, 
                 self.unmapped_inds))
    
    def get_padding(self):
        if len(self.mapped_inds) == 0:
            return 0, 0
        if len(self.unmapped_inds) == 0:
            return 0, 0
        # pad left
        pad_left = 0
        for i in xrange(self.start_ind, self.mapped_inds[0]):
            pad_left += len(self.split_dict[i][0].seq)
        # pad right
        pad_right = 0
        for i in xrange(self.mapped_inds[-1] + 1, self.end_ind):
            pad_right += len(self.split_dict[i][0].seq)
        return (pad_left, pad_right) if self.strand == 0 else (pad_right, pad_left)
    

def find_split_points(refmaps, num_clusters):
    cluster_intervals = []
    # build a mapping between splits in the refmaps and the clusters 
    # those splits are associated with.  we call this an 'ind_clust_map'
    # which associates a split index with a set of intervals
    # that represent clusters    
    ind_clust_map = collections.defaultdict(lambda: set())
    for refmap in refmaps.itervalues():
        for interval in refmap.get_clusters():
            # add to master list of read clusters
            clust_id = len(cluster_intervals)
            cluster_intervals.append(interval)
            # build an index from cluster index to id of reference/strand
            # object containing the ReadClusters at that index
            ind_rclust_dict = interval.value 
            for ind in ind_rclust_dict:
                ind_clust_map[ind].add(clust_id)    
    # now walk through cluster indexes in order to find split points
    start_ind = 0
    mapped_inds = []
    current_clust_ids = set()
    concordant_clust_ids = []
    for split_ind in xrange(num_clusters):
        # if there are no clusters then this index
        # is unmapped and not useful
        if split_ind not in ind_clust_map:
            continue  
        # get items at this index
        clust_ids = set(ind_clust_map[split_ind])
        if len(current_clust_ids) == 0:
            # initialize the current clusters to the 
            # first mapped index
            current_clust_ids.update(clust_ids)
        elif current_clust_ids.isdisjoint(clust_ids):
            # save the set of concordant clusters and the split index
            concordant_clust_ids.append((start_ind, split_ind, mapped_inds, current_clust_ids))
            # initialize new clusters
            current_clust_ids = clust_ids            
            start_ind = mapped_inds[-1] + 1
            mapped_inds = []
        else:
            # intersect index clusters together to 
            # reduce number of current clusters
            current_clust_ids.intersection_update(clust_ids)
        # keep track of the last mapped index so that any
        # unmapped indexes in between can be attributed to 
        # both the 5' and 3' split genes
        mapped_inds.append(split_ind)
    # add the last set of clusters
    if len(current_clust_ids) > 0:
        concordant_clust_ids.append((start_ind, split_ind+1, mapped_inds, current_clust_ids))
    return cluster_intervals, concordant_clust_ids

def find_split_read_clusters(partition_splits, max_indel_size, max_multihits):
    refmaps = {}
    split_mapping_codes = []
    unmapped_read_dict = collections.defaultdict(lambda: [])
    for split_ind, split_reads in enumerate(partition_splits):
        mapping_codes = set()                         
        for r in split_reads:
            #print 'IND', split_ind, 'mapping code', _mapping_code_strings[get_mapping_code(r, max_multihits)]            
            # keep track of mapping results for reads in this split
            mapping_codes.add(get_mapping_code(r, max_multihits))
            if r.is_unmapped:
                unmapped_read_dict[split_ind].append(r)
                continue
            # cluster reads by reference name and strand
            strand = int(r.is_reverse)            
            rkey = (r.rname, strand)
            if rkey not in refmaps:
                refmaps[rkey] = RefCluster(r.rname, strand, max_indel_size)
            refmaps[rkey].add(r.pos, r.aend, split_ind, r)
        split_mapping_codes.append(mapping_codes)
    # convert unmapped reads index dict to a static dict
    unmapped_read_dict = dict(unmapped_read_dict)
    # search reference maps by index to find split points
    cluster_intervals, concordant_clust_ids = find_split_points(refmaps, len(partition_splits))
    # convert from cluster ids to read cluster intervals
    concordant_clusters = []
    for start_ind, end_ind, mapped_inds, clust_ids in concordant_clust_ids:
        print 'SR START', start_ind, 'END', end_ind, 'MAPPED INDS', mapped_inds, 'IDS', clust_ids
        # determine indexes of unmapped splits in this cluster and 
        # get lists of unmapped reads        
        unmapped_inds = set(xrange(start_ind, end_ind)).difference(mapped_inds)
        rclusts = []        
        for id in clust_ids:
            interval = cluster_intervals[id]
            # build a new ReadCluster object with complete information about 
            # start,end indexes and which indexes are mapping/nonmapping.  
            # thus clusters will now be contiguous lists of split indexes with 
            # padding information at the start/end
            rclust = ReadCluster(interval, start_ind, end_ind, mapped_inds, 
                                 unmapped_inds, unmapped_read_dict)
            rclusts.append(rclust)
        concordant_clusters.append((start_ind, end_ind, mapped_inds, unmapped_inds, tuple(rclusts)))
#        if len(partition_splits) < 3:
#            continue
#        for id in clust_ids:
#            c = rclusters[id]
#            print 'CLUSTER'
#            for i in xrange(start_ind, end_ind):
#                for r in c.split_read_dict[i]:
#                    print r
    # fill in unmapped reads within cluster
    #rclust.update_missing(unmapped_read_dict)        
    #print 'UNMAPPED', dict(unmapped_read_dict)
    #print 'IND', split_ind, split_mapping_codes[split_ind], ind_clust_map[split_ind]
    return split_mapping_codes, concordant_clusters

def merge_read_clusters(clusters1, clusters2, max_dist, library_type):
    assert library_type[0] != library_type[1]
    # TODO: currently only support "fr" or "rf" libraries
    paired_clusters = list(itertools.chain(clusters1, reversed(clusters2)))
    refmaps = {}
    for cluster_ind, cluster_info in enumerate(paired_clusters):
        # figure out which mate this is (0 if index is less than
        # number of clusters in read1, 1 otherwise)
        mate = cluster_ind >= len(clusters1)
        # add the ReadCluster objects to cluster trees
        start_ind, end_ind, mapped_inds, unmapped_inds, rclusts = cluster_info            
        for rclust in rclusts:
            # annotate each cluster with mate
            rclust.mate = mate
            # we reverse strand of one of the mates, and here arbitrarily
            # choose the second mate
            strand = int(not rclust.strand) if mate == 1 else rclust.strand
            # use the flipped strand in the key - we must dig into the
            # individual ReadClusters to get the original strands back
            rkey = (rclust.rname, strand)
            if rkey not in refmaps:
                refmaps[rkey] = RefCluster(rclust.rname, strand, max_dist)
            refmaps[rkey].add(rclust.start, rclust.end, cluster_ind, rclust)
    # organize clusters based on the indexes they incorporate
    ind_cluster_dict = collections.defaultdict(lambda: [])
    for rkey, refmap in refmaps.iteritems():
        for interval in refmap.get_clusters():
            # build an index from cluster indexes to cluster interval data
            ind_dict = interval.value
            ind_cluster_dict[tuple(sorted(ind_dict))].append(interval)
    # sorted index tuples from low to hi to ensure reading from read1 -> read2
    # if read1 and read2 overlap, this should still work.
    paired_clusters = []
    for ind_tuple in sorted(ind_cluster_dict):
        print 'TUPLE', ind_tuple
        paired_clusters.append(ind_cluster_dict[ind_tuple])
    return paired_clusters


class DiscordantCluster(object):
    __slots__ = ('rname', 'start', 'end', 'strand', 'pad_start', 'pad_end', 
                 'multimaps', 'seq', 'qual')
    def __init__(self, rname, start, end, strand, pad_start, pad_end, multimaps,
                 seq=None, qual=None):
        self.rname = rname
        self.start = start
        self.end = end
        self.strand = strand
        self.pad_start = pad_start
        self.pad_end = pad_end
        self.multimaps = multimaps
        self.seq = seq
        self.qual = qual
    def __repr__(self):
        return ("<%s(rname=%s, strand=%s, start=%d, end=%d, pad_start=%d, pad_end=%d, multimaps=%d, seq=%s, qual=%s)>" %
                (self.__class__.__name__, self.rname, self.strand, self.start, self.end,
                 self.pad_start, self.pad_end, self.multimaps, self.seq, self.qual))
    def to_list(self):
        return [self.rname, self.start, self.end, self.strand, 
                self.pad_start, self.pad_end, self.multimaps, 
                self.seq, self.qual]
    @staticmethod
    def from_line(self, line):
    

class DiscordantFragment(object):
    __slots__ = ('qname', 'discordant_type', 'read1_is_5prime', 'clust1', 'clust2')
    NONMAPPING = 0
    CONCORDANT = 1
    DISCORDANT_INNER = 2
    DISCORDANT_READ1 = 3
    DISCORDANT_READ2 = 4
    DISCORDANT_OVERLAPPING = 5
    DISCORDANT_GENOME = 6
    DISCORDANT_COMPLEX = 7
    _discordant_types = ["NONMAPPING",
                         "CONCORDANT",                         
                         "DISCORDANT_INNER",
                         "DISCORDANT_READ1",
                         "DISCORDANT_READ2",
                         "DISCORDANT_OVERLAPPING",
                         "DISCORDANT_GENOME",
                         "DISCORDANT_COMPLEX"]

    @staticmethod
    def get_discordant_type(clusters1, clusters2, paired_clusters):
        if len(paired_clusters) == 0:
            discordant_type = DiscordantFragment.NONMAPPING
        elif len(paired_clusters) == 1:
            discordant_type = DiscordantFragment.CONCORDANT
        elif len(paired_clusters) > 2:
            discordant_type = DiscordantFragment.DISCORDANT_COMPLEX            
        elif len(paired_clusters) == 2:
            if len(clusters1) > 1 and len(clusters2) > 1:
                discordant_type = DiscordantFragment.DISCORDANT_OVERLAPPING
            elif len(clusters1) > 1:
                discordant_type = DiscordantFragment.DISCORDANT_READ1
            elif len(clusters2) > 1:
                discordant_type = DiscordantFragment.DISCORDANT_READ2
            else:
                discordant_type = DiscordantFragment.DISCORDANT_INNER
        return discordant_type
        
    def __init__(self, qname, discordant_type, read1_is_5prime,
                 clust1, clust2): 
        self.qname = qname
        self.discordant_type = discordant_type
        self.read1_is_5prime = read1_is_5prime
        self.clust1 = clust1
        self.clust2 = clust2

    def to_list(self):
        return [self.qname, self._discordant_types[self.discordant_type], 
                self.read1_is_5prime] + self.clust1.to_list() + \
                self.clust2.to_list()

#def get_padded_bounds(interval):
#    pad_start, pad_end = interval.start, interval.end
#    ind_rclust_dict = interval.value        
#    for rclusts in ind_rclust_dict.itervalues():
#        for rclust in rclusts:
#            left, right = rclust.get_padding()
#            #print 'PAD', pad_start, pad_end, 'RCLUST', rclust, 'LEFT, RIGHT', left, right
#            if (rclust.start - left) < pad_start:
#                pad_start = rclust.start - left
#            if (rclust.end + right) > pad_end:
#                pad_end = rclust.end + right
#    return pad_start, pad_end

def interval_to_discordant_cluster(interval, tid_rname_map):
    chrom = "*" if interval.chrom == -1 else tid_rname_map[interval.chrom]    
    strand = "+" if interval.strand == 0 else "-"
    pad_start, pad_end = interval.start, interval.end
    multimaps = None
    ind_rclust_dict = interval.value        
    for rclusts in ind_rclust_dict.itervalues():
        for rclust in rclusts:
            left, right = rclust.get_padding()
            #print 'PAD', pad_start, pad_end, 'RCLUST', rclust, 'LEFT, RIGHT', left, right
            if (rclust.start - left) < pad_start:
                pad_start = rclust.start - left
            if (rclust.end + right) > pad_end:
                pad_end = rclust.end + right
            if (multimaps is None) or multimaps < rclust.multimaps:
                multimaps = rclust.multimaps
    clust = DiscordantCluster(chrom, 
                              interval.start, 
                              interval.end, 
                              strand,
                              pad_start, 
                              pad_end, 
                              multimaps)
    return clust


def clusters_to_discordant_pairs(qname, clusters1, clusters2, paired_clusters, gene_genome_map, 
                                 tid_rname_map):
    # if there are less than 2 clusters, then the reads are not discordant
    # according to the constraints set in the program
    discordant_type = DiscordantFragment.get_discordant_type(clusters1, 
                                                             clusters2, 
                                                             paired_clusters)
    print 'RESULT', DiscordantFragment._discordant_types[discordant_type], len(paired_clusters)
    if len(paired_clusters) != 2:
        return []
    # find candidate 5' and 3' clusters. for genes organize as 
    # 5' and 3'.  for genome hits keep as read1 and read2
    gene_pairs = (([],[]),([],[]))
    genome_pairs = ([],[])    
    for clust_num, intervals in enumerate(paired_clusters):
        # first item in paired clusters will guarantee to contain
        # read1 and strand == 0 means read1 is 5' (strand == 1 means
        # read1 is 3'
        for interval in intervals:                
            print 'INTERVAL', clust_num, 'STRAND', interval.strand, interval
            # convert interval reference ID to chromosome name
            if gene_genome_map[interval.chrom] is None:
                genome_pairs[clust_num].append(interval)
            else:
                mate = clust_num if interval.strand == 0 else int(not clust_num)
                gene_pairs[interval.strand][mate].append(interval)
    #print 'GENE PAIRS strand=0', gene_pairs[0]
    #print 'GENE PAIRS strand=1', gene_pairs[1]
    discordant_pairs = []
    # join 5'/3' gene hits
    for strand, mate_intervals in enumerate(gene_pairs):
        intervals5p, intervals3p = mate_intervals
        for interval5p in intervals5p:
            clust5p = interval_to_discordant_cluster(interval5p, tid_rname_map)            
            for interval3p in intervals3p:
                clust3p = interval_to_discordant_cluster(interval3p, tid_rname_map)            
                # combine 5'/3' partners
                discordant_pairs.append(DiscordantFragment(qname, discordant_type, (strand == 0), clust5p, clust3p))
    # if not gene hits, then output genome hits
    if len(discordant_pairs) == 0:
        intervals1 = genome_pairs[0]
        intervals2 = genome_pairs[1]
        for interval1 in intervals1:
            clust1 = interval_to_discordant_cluster(interval1, tid_rname_map)            
            for interval2 in intervals2:
                clust2 = interval_to_discordant_cluster(interval2, tid_rname_map)            
                # combine read1/read2 intervals
                discordant_type = DiscordantFragment.DISCORDANT_GENOME
                discordant_pairs.append(DiscordantFragment(qname, discordant_type, True, clust1, clust2))
    return discordant_pairs


def find_discordant_pairs(pe_reads, tid_rname_map, gene_genome_map, max_indel_size, 
                          max_isize, max_multihits, library_type):
    pe_clusters = ([],[])
    all_mapping_codes = set()
    qname = pe_reads[0][0][0][0].qname    
    # search for discordant read clusters in individual reads first
    # before using paired-end information
    for mate, partitions in enumerate(pe_reads):
        for partition_ind, partition_splits in enumerate(partitions):
            # check for discordant reads 
            split_mapping_codes, read_clusters = find_split_read_clusters(partition_splits, max_indel_size, max_multihits)
            # update aggregated mapping codes
            all_mapping_codes.update(*split_mapping_codes)
            # store mate clusters
            pe_clusters[mate].append(read_clusters)
    # quit early if both mates unmapped
    if (len(all_mapping_codes) == 1 and MAPPING not in all_mapping_codes):
        return
    # combine paired-end cluster information to predict discordant reads
    pairs = []
    for read1_clusters in pe_clusters[0]:
        for read2_clusters in pe_clusters[1]:
            # try to combine 5'/3' partners
            paired_clusters = merge_read_clusters(read1_clusters, 
                                                  read2_clusters, 
                                                  max_isize, 
                                                  library_type)
            # make discordant pair objects
            pairs.extend(clusters_to_discordant_pairs(qname, 
                                                      read1_clusters, 
                                                      read2_clusters, 
                                                      paired_clusters,                                                       
                                                      gene_genome_map,
                                                      tid_rname_map))
    return pairs

def find_discordant_reads(bamfh, output_file, 
                          gene_genome_map, max_indel_size, 
                          max_isize, max_multihits, library_type):
    refs = bamfh.references
    outfh = open(output_file, "w")
    for pe_reads in parse_reads(bamfh):
        pairs = find_discordant_pairs(pe_reads, refs, gene_genome_map, 
                                      max_indel_size, max_isize, 
                                      max_multihits, library_type)
        for pair in pairs:
            print >>outfh, '\t'.join(map(str, pair.to_list()))
    outfh.close()
        
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
    output_file = args[1]
    gene_file = os.path.join(options.index_dir, config.GENE_FEATURE_FILE)
    library_type = parse_library_type(options.library_type)
    # open bam file
    bamfh = pysam.Samfile(input_bam_file, "rb")
    # build genome map
    logging.info("Loading gene table")
    gene_genome_map, gene_trees = build_gene_maps(bamfh, gene_file)    
    find_discordant_reads(bamfh, output_file, gene_genome_map,
                          max_indel_size=options.max_indel_size, 
                          max_isize=options.max_fragment_length,
                          max_multihits=options.multihits,
                          library_type=library_type)
    bamfh.close()

if __name__ == '__main__':
    main()
