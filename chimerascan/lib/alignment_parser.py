'''
Created on Jan 22, 2011

@author: mkiyer
'''
from base import SamTags

def parse_sr_sam_file(samfh):    
    reads = []
    for read in samfh:        
        if len(reads) > 0 and read.qname != reads[-1].qname:
            yield reads
            reads = []
        reads.append(read)
    if len(reads) > 0:
        yield reads

def parse_segmented_pe_sam_file(bamfh):
    pe_reads = ([], [])
    # reads must be binned by qname, mate, hit, and segment
    # so initialize to mate 0, hit 0, segment 0
    num_reads = 0
    prev_qname = None
    for read in bamfh:
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

def parse_segmented_sr_sam_file(bamfh):
    reads = []
    # reads must be binned by qname, mate, hit, and segment
    # so initialize to mate 0, hit 0, segment 0
    num_reads = 0
    prev_qname = None
    for read in bamfh:
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
            yield reads
            # reset state variables
            reads = []
            num_reads = 0
        prev_qname = qname
        # initialize hits
        if len(reads) == 0:
            reads.extend([list() for x in xrange(num_split_partitions)])
        # initialize hit segments
        if len(reads[partition_ind]) == 0:
            reads[partition_ind].extend([list() for x in xrange(num_splits)])
        split_reads = reads[partition_ind][split_ind]
        # initialize segment mappings
        if len(split_reads) == 0:
            split_reads.extend([None for x in xrange(num_mappings)])
        # add segment to hit/mate/read
        split_reads[mapping_ind] = read
        num_reads += 1
    if num_reads > 0:
        yield reads
