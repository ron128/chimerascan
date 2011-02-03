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
from chimerascan import pysam
from chimerascan.bx.intersection import Interval, IntervalTree

# local imports
from chimerascan.lib import config
from chimerascan.lib.base import parse_library_type, parse_bool, parse_string_none, SamTags
from chimerascan.lib.gene_to_genome import build_gene_maps

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


class RefCluster(object):
    def __init__(self, rname, strand, max_dist):
        self.rname = rname
        self.strand = strand
        self.vals = []
        self.max_dist = max_dist

    def add(self, start, end, ind, val):
        self.vals.append((start, end, ind, val))

    def get_clusters(self):
        if len(self.vals) == 0:
            return
        # sort values in order of increasing start position
        self.vals.sort(key=operator.itemgetter(0))
        # initialize cluster
        valiter = iter(self.vals)
        start, end, ind, val = valiter.next()
        ind_val_dict = collections.defaultdict(lambda: [])
        ind_val_dict[ind].append(val)
        # find clusters
        for next_start, next_end, ind, val in valiter:
            if next_start > (end + self.max_dist):
                # this interval is outside bounds of current cluster
                yield Interval(start, end, chrom=self.rname, strand=self.strand, value=ind_val_dict)
                ind_val_dict = collections.defaultdict(lambda: [])
                start = next_start
            # update values
            if next_end > end:
                end = next_end
            ind_val_dict[ind].append(val)
        if len(ind_val_dict) > 0:
            yield Interval(start, end, chrom=self.rname, strand=self.strand, value=ind_val_dict)

class ReadCluster(object):
    def __init__(self, interval, start_ind, end_ind, mapped_inds, 
                 unmapped_inds, unmapped_read_dict):
        # access interval object with start, end, strand, chrom info
        self.start = interval.start
        self.end = interval.end
        self.rname = interval.chrom
        self.strand = interval.strand
        self.split_dict = {}
        self.multimaps = None
        # copy split dict from clustered intervals (value)
        for ind,reads in interval.value.iteritems():
            self.split_dict[ind] = reads
            # find minimum multimaps for all reads in cluster
            min_multimaps = min(r.opt('NH') for r in reads)
            if (self.multimaps is None) or self.multimaps < min_multimaps:
                self.multimaps = min_multimaps
        # set to zero if still None
        if self.multimaps is None:
            self.multimaps = 0
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
        # set mate to 0 (read1) or 1 (read2)
        self.mate = self.split_dict[start_ind][0].is_read2

    def __repr__(self):
        return ("<%s(rname=%s, strand=%s, start=%d, end=%d, multimaps=%d, mate=%d, start_ind=%d, "
                "end_ind=%d, mapped_inds=%s, unmapped_inds=%s)>" % 
                (self.__class__.__name__, self.rname, self.strand, self.start, self.end,
                 self.multimaps, self.mate, self.start_ind, self.end_ind, self.mapped_inds, 
                 self.unmapped_inds))
    
    @property
    def is_unmapped(self):
        return self.rname == -1
    
    def get_padding(self, padding=0):
        '''
        padding: extra padding to add.  used when reads were trimmed during
        mapping and want to see whether they could be spanning
        '''
        if self.strand == 0:
            pad_left, pad_right = 0, padding
        else:
            pad_left, pad_right = padding, 0
        if (len(self.mapped_inds) == 0 or
            len(self.unmapped_inds) == 0):
            return pad_left, pad_right
        # pad left
        for i in xrange(self.start_ind, self.mapped_inds[0]):
            pad_left += len(self.split_dict[i][0].seq)
        # pad right
        for i in xrange(self.mapped_inds[-1] + 1, self.end_ind):
            pad_right += len(self.split_dict[i][0].seq)
        # add extra padding
        if self.strand == 0:
            return (pad_left, pad_right)
        else:
            return (pad_right, pad_left)


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
                #print 'UNMAPPED', 'IND', split_ind, 'READ', r
                unmapped_read_dict[split_ind].append(r)
                continue
            #print 'MAPPED', 'IND', split_ind, 'READ', r
            # cluster reads by reference name and strand
            strand = int(r.is_reverse)            
            rkey = (r.rname, strand)
            if rkey not in refmaps:
                refmaps[rkey] = RefCluster(r.rname, strand, max_indel_size)
            refmaps[rkey].add(r.pos, r.aend, split_ind, r)
        split_mapping_codes.append(mapping_codes)
    # convert unmapped reads index dict to a static dict
    unmapped_read_dict = dict(unmapped_read_dict)
    if all((MAPPING not in codes) for codes in split_mapping_codes):
        # if there are no mapped reads, create a dummy "unmapped" ReadCluster
        # and return early
        interval = Interval(0, 0, chrom=-1, strand=0, value={})        
        rclust = ReadCluster(interval, start_ind=0, end_ind=len(partition_splits), 
                             mapped_inds=[], unmapped_inds=set(unmapped_read_dict),
                             unmapped_read_dict=unmapped_read_dict)
        # results are returned as a list of tuples, where each tuple represents
        # a set of reads that cluster together. return a singleton tuple here  
        return split_mapping_codes, [(rclust,)]
    # search reference maps by index to find split points
    cluster_intervals, concordant_clust_ids = find_split_points(refmaps, len(partition_splits))
    # convert from cluster ids to read cluster intervals
    concordant_clusters = []
    for start_ind, end_ind, mapped_inds, clust_ids in concordant_clust_ids:
        #print 'SR START', start_ind, 'END', end_ind, 'MAPPED INDS', mapped_inds, 'IDS', clust_ids
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
        concordant_clusters.append(tuple(rclusts))
    return split_mapping_codes, concordant_clusters

def merge_read_clusters(clusters1, clusters2, max_dist, library_type):
    assert library_type[0] != library_type[1]
    # TODO: currently only support "fr" or "rf" libraries
    paired_clusters = list(itertools.chain(clusters1, reversed(clusters2)))
    refmaps = {}    
    for cluster_ind, rclusts in enumerate(paired_clusters):
        # add the ReadCluster objects to cluster trees
        for rclust in rclusts:
            # ignore unmapped read clusters
            if rclust.is_unmapped:
                continue
            # we reverse strand of one of the mates, and by convention
            # choose the second mate
            strand = int(not rclust.strand) if rclust.mate == 1 else rclust.strand
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
    def cmplen(x,y):
        res = -1*cmp(len(x),len(y))
        if res != 0: return res
        return cmp(x,y)
    paired_clusters = []
    used_inds = set()
    for ind_tuple in sorted(ind_cluster_dict, cmp=cmplen):
        if used_inds.isdisjoint(ind_tuple):            
            used_inds.update(ind_tuple)
            paired_clusters.append(ind_cluster_dict[ind_tuple])
        # TODO: analyze complex mapping tuples further
        #else:
        #    logging.warning("TUPLE %s NOT DISJOINT (ALL TUPLES=%s" % (ind_tuple, sorted(ind_cluster_dict, cmp=cmplen)))
    return paired_clusters


class DiscordantCluster(object):
    __slots__ = ('rname', 'start', 'end', 'strand', 'pad_start', 'pad_end', 
                 'multimaps', 'seq', 'qual')
    def __init__(self, rname, start, end, strand, pad_start, pad_end, multimaps,
                 seq=None, 
                 qual=None):
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
        return ("<%s(rname=%s, strand=%s, start=%d, end=%d, pad_start=%d, "
                "pad_end=%d, multimaps=%d, seq=%s, qual=%s)>" %
                (self.__class__.__name__, self.rname, self.strand, self.start, self.end,
                 self.pad_start, self.pad_end, self.multimaps, self.seq, self.qual))
    def to_list(self):
        return [self.rname, self.start, self.end, self.strand, 
                self.pad_start, self.pad_end, self.multimaps,
                self.seq, self.qual]
    @staticmethod
    def from_list(fields):        
        seq = parse_string_none(fields[7])
        qual = parse_string_none(fields[8])
        return DiscordantCluster(fields[0], int(fields[1]), int(fields[2]), 
                                 fields[3], int(fields[4]), int(fields[5]),
                                 int(fields[6]), seq, qual)

    
class DiscordantFragment(object):
    __slots__ = ('qname', 'is_genome', 'discordant5p', 'discordant3p',
                 'code', 'read1_is_sense', 'clust5p', 'clust3p')
    # columns in the tabular output that contain the references
    # (useful for sorting)
    REF1_COL = 6
    REF2_COL = 15
    
    # fragment type codes
    NA = 0
    NONMAPPING = 1
    CONCORDANT_SINGLE = 2
    DISCORDANT_SINGLE = 3
    DISCORDANT_SINGLE_COMPLEX = 4
    CONCORDANT_PAIRED = 5
    DISCORDANT_PAIRED = 6
    DISCORDANT_PAIRED_COMPLEX = 7
    _discordant_codes = ["NA", 
                         "NONMAPPING",
                         "CONCORDANT_SINGLE",
                         "DISCORDANT_SINGLE",
                         "DISCORDANT_SINGLE_COMPLEX",
                         "CONCORDANT_PAIRED",
                         "DISCORDANT_PAIRED",
                         "DISCORDANT_PAIRED_COMPLEX"]
    # flag bits
    FLAG_DISCORDANT_5P = 0b10
    FLAG_DISCORDANT_3P = 0b1

    def __init__(self, qname, read1_is_sense, clust5p, clust3p,
                 code=0, is_genome=False, discordant5p=False, 
                 discordant3p=False): 
        self.qname = qname
        self.read1_is_sense = read1_is_sense
        self.clust5p = clust5p
        self.clust3p = clust3p
        self.code = code
        self.is_genome = is_genome
        self.discordant5p = discordant5p
        self.discordant3p = discordant3p

    def __repr__(self):
        return ("<%s(qname=%s, read1_is_sense=%s, clust5p=%s, clust3p=%s "
                "is_genome=%s, discordant5p=%s, discordant3p=%s, "
                "code=%d, string_code=%s)>" %
                (self.__class__.__name__, self.qname, self.read1_is_sense,
                 self.clust5p, self.clust3p, 
                 self.is_genome, self.discordant5p, self.discordant3p, 
                 self.code, self._discordant_codes[self.code]))

    @property
    def clust1(self):
        return self.clust5p if self.read1_is_sense else self.clust3p
    @property
    def clust2(self):
        return self.clust3p if self.read1_is_sense else self.clust5p

    def set_flags(self, nclusts1, nclusts2, nclustspe):
        # set 5'/3' discordant flags according to sense/antisense
        # orientation and cluster mappings
        if self.read1_is_sense:
            self.discordant5p = (nclusts1 > 1)
            self.discordant3p = (nclusts2 > 1)
        else:
            self.discordant5p = (nclusts2 > 1)
            self.discordant3p = (nclusts1 > 1)                    
        if nclustspe == 0:
            self.code = self.NONMAPPING
        elif nclustspe == 1:
            if (nclusts1 == 0) or (nclusts2 == 0):
                # one of the reads is concordant and the other is unmapped
                self.code = self.CONCORDANT_SINGLE
            else:
                # both reads are mapped and concordant
                self.code = self.CONCORDANT_PAIRED
        elif nclustspe == 2:
            if (nclusts1 > 1) or (nclusts2 > 1):
                # one of the reads is discordant and the other is
                # unmapped
                if (nclusts1 > 2) or (nclusts2 > 2):        
                    self.code = self.DISCORDANT_SINGLE_COMPLEX
                else:
                    self.code = self.DISCORDANT_SINGLE
            else:
                self.code = self.DISCORDANT_PAIRED
        else:
            self.code = self.DISCORDANT_PAIRED_COMPLEX
        return self

    def _pack_flags(self):
        '''pack flags into a single integer value for reading/writing to file'''
        flags = ((self.discordant5p << self.FLAG_DISCORDANT_5P) |
                 (self.discordant3p << self.FLAG_DISCORDANT_3P))
        return flags
    
    def _unpack_flags(self, val):
        '''unpack flags and set attributes'''
        self.discordant5p = bool(val & self.FLAG_DISCORDANT_5P)
        self.discordant3p = bool(val & self.FLAG_DISCORDANT_3P)

    def to_list(self):
        genome = "GENOME" if self.is_genome else "GENE"
        return ([self.qname, int(self.read1_is_sense), genome,
                 self._discordant_codes[self.code], self._pack_flags()] +
                self.clust5p.to_list() + self.clust3p.to_list())

    @staticmethod
    def from_list(fields):
        qname = fields[0]        
        read1_is_sense = bool(int(fields[1]))
        is_genome = True if fields[2] == "GENOME" else False
        code = DiscordantFragment._discordant_codes.index(fields[3])
        flags = int(fields[4])
        clust5p = DiscordantCluster.from_list(fields[5:14])
        clust3p = DiscordantCluster.from_list(fields[14:23])
        f = DiscordantFragment(qname, read1_is_sense, clust5p, clust3p,
                               code=code, is_genome=is_genome)
        f._unpack_flags(flags)
        return f
        

def interval_to_discordant_cluster(interval, tid_rname_map, gene_genome_map,
                                   padding):
    if interval.chrom == -1:
        chrom = "*"
    elif gene_genome_map[interval.chrom] is not None:
        chrom = gene_genome_map[interval.chrom].tx_name
    else:
        chrom = tid_rname_map[interval.chrom]    
    strand = "+" if interval.strand == 0 else "-"
    pad_start, pad_end = interval.start, interval.end
    multimaps = None
    ind_rclust_dict = interval.value    
    for rclusts in ind_rclust_dict.itervalues():
        for rclust in rclusts:
            left, right = rclust.get_padding(padding=padding)
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

def count_cluster_mappings(clusts):
    '''
    return number of mapped ReadCluster objects within the
    'clusts' list
    '''
    num_mapped = 0
    for rclusts in clusts:
        is_mapped = 0
        for rclust in rclusts:
            if rclust.is_unmapped:
                continue
            is_mapped = 1
            break
        num_mapped += is_mapped
    return num_mapped

def clusters_to_discordant_fragments(qname, clusters1, clusters2, 
                                     paired_clusters, gene_genome_map, 
                                     tid_rname_map, padding):
    assert len(paired_clusters) > 1
    nclusts1 = count_cluster_mappings(clusters1)
    nclusts2 = count_cluster_mappings(clusters2)
    nclustspe = len(paired_clusters)        
    # bin clusters as genes/genome.  enforce strandedness for 
    # gene clusters  
    gene_clusts = (([], []), ([], []))
    genome_clusts = ([],[])
    for clust_num, intervals in enumerate(paired_clusters):
        for interval in intervals:
            # create DiscordantCluster objects and bin by cluster index and
            # strand (genes only) 
            discordant_clust = \
                interval_to_discordant_cluster(interval, tid_rname_map, 
                                               gene_genome_map, padding)
            if gene_genome_map[interval.chrom] is None:
                genome_clusts[clust_num].append(discordant_clust)
            else:
                gene_clusts[clust_num][interval.strand].append(discordant_clust)
    # paired cluster hits
    frags = []
    for strand in (0,1):
        # find correct strand for discordant clusters that have
        # read1 or read2 unmapped
        if nclusts1 == 0:
            read1_is_sense = (strand == 1)
        else:
            read1_is_sense = (strand == 0)        
        clusts1, clusts2 = gene_clusts[0][strand], gene_clusts[1][strand]
        for clust1 in clusts1:
            for clust2 in clusts2:
                if read1_is_sense:
                    clust5p, clust3p = clust1, clust2
                else:
                    clust5p, clust3p = clust2, clust1                
                f = DiscordantFragment(qname, read1_is_sense, clust5p, clust3p,
                                       is_genome=False)
                f.set_flags(nclusts1, nclusts2, nclustspe)
                frags.append(f)                
    if len(frags) > 0:
        return frags, []
    # if no gene hits, then find genome hits
    # correct for strand if read1 is unmapped
    clusts1, clusts2 = genome_clusts[0], genome_clusts[1]
    for clust1 in clusts1:
        for clust2 in clusts2:
            if nclusts1 == 0:
                read1_is_sense = int(clust2.strand == 1)
            else:
                read1_is_sense = int(clust1.strand == 0)
            f = DiscordantFragment(qname, read1_is_sense, clust1, clust2,
                                   is_genome=False)
            f.set_flags(nclusts1, nclusts2, nclustspe)
            frags.append(f)                
    return [], frags


def find_discordant_pairs(pe_reads, tid_rname_map, gene_genome_map, 
                          max_indel_size, max_isize, max_multihits, 
                          library_type, padding):
    pe_clusters = ([],[])
    mate_mapping_codes = (set(), set())
    qname = pe_reads[0][0][0][0].qname
    #print 'QNAME', qname 
    # search for discordant read clusters in individual reads first
    # before using paired-end information
    for mate, partitions in enumerate(pe_reads):
        #print 'MATE', mate, 'PARTITIONS', len(partitions)
        for partition_ind, partition_splits in enumerate(partitions):
            #print 'PARTITION', partition_ind, 'SPLITS', len(partition_splits)
            # check for discordant reads 
            split_mapping_codes, read_clusters = \
                find_split_read_clusters(partition_splits, max_indel_size, 
                                         max_multihits)
            # update aggregated mapping codes
            mate_mapping_codes[mate].update(*split_mapping_codes)
            # store mate clusters
            pe_clusters[mate].append(read_clusters)
    # find out how many clusters we produced from each read in the pair
    #print 'QNAME', qname, 'MAPPING CODES', mate_mapping_codes
    both_unmapped = ((MAPPING not in mate_mapping_codes[0]) and
                     (MAPPING not in mate_mapping_codes[1]))
    # if unmapped return early
    if both_unmapped:
        return [],[]
    # combine paired-end cluster information to predict discordant reads
    gene_frags = []
    genome_frags = []
    for read1_clusters in pe_clusters[0]:
        for read2_clusters in pe_clusters[1]:
            # combine 5'/3' partners
            paired_clusters = merge_read_clusters(read1_clusters, 
                                                  read2_clusters, 
                                                  max_isize, 
                                                  library_type)
            if len(paired_clusters) <= 1:
                # if the reads cluster together than they are
                # not discordant
                continue
            if len(paired_clusters) > 2:
                # if the reads form more than two clusters they
                # are "complex"
                continue
            # make discordant pair objects
            gene_hits, genome_hits = \
                clusters_to_discordant_fragments(qname, 
                                                 read1_clusters, 
                                                 read2_clusters, 
                                                 paired_clusters,                                                       
                                                 gene_genome_map,
                                                 tid_rname_map,
                                                 padding)
            gene_frags.extend(gene_hits)
            genome_frags.extend(genome_hits)
    return gene_frags, genome_frags


def find_discordant_reads(bamfh, gene_output_file, genome_output_file, 
                          gene_file, max_indel_size, max_isize, 
                          max_multihits, library_type, padding):
    # build genome map
    logging.info("Loading gene table")
    gene_genome_map, gene_trees = build_gene_maps(bamfh, gene_file)    
    logging.info("Finding discordant reads")
    refs = bamfh.references
    gene_fh = open(gene_output_file, "w")
    genome_fh = open(genome_output_file, "w")
    for pe_reads in parse_reads(bamfh):    
        gene_pairs, genome_pairs = \
            find_discordant_pairs(pe_reads, refs, gene_genome_map, 
                                  max_indel_size, max_isize, 
                                  max_multihits, library_type,
                                  padding)
        for pair in gene_pairs:
            print >>gene_fh, '\t'.join(map(str, pair.to_list()))
        for pair in genome_pairs:
            print >>genome_fh, '\t'.join(map(str, pair.to_list()))
    gene_fh.close()
    genome_fh.close()
        
def main():
    from optparse import OptionParser
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = OptionParser("usage: %prog [options] <bam> <out.bedpe>")
    parser.add_option("--index", dest="index_dir",
                      help="Path to chimerascan index directory")
    parser.add_option('--max-fragment-length', dest="max_fragment_length", 
                      type="int", default=1000)
    parser.add_option('--max-indel-size', dest="max_indel_size", 
                      type="int", default=100)
    parser.add_option('--library-type', dest="library_type", default="fr")
    parser.add_option('--multihits', type="int", default=1)
    parser.add_option('--padding', type="int", default=0)
    options, args = parser.parse_args()
    input_bam_file = args[0]
    gene_output_file = args[1]
    genome_output_file = args[2]    
    gene_feature_file = os.path.join(options.index_dir, config.GENE_FEATURE_FILE)
    library_type = parse_library_type(options.library_type)
    # open bam file
    bamfh = pysam.Samfile(input_bam_file, "rb")
    find_discordant_reads(bamfh, gene_output_file, genome_output_file,                           
                          gene_feature_file,
                          max_indel_size=options.max_indel_size, 
                          max_isize=options.max_fragment_length,
                          max_multihits=options.multihits,
                          library_type=library_type,
                          padding=options.padding)
    bamfh.close()

if __name__ == '__main__':
    main()
