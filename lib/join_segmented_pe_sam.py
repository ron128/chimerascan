'''
Created on Jan 7, 2011

@author: mkiyer
'''
import sys
import tempfile
import os
import logging
import collections
import subprocess
import multiprocessing
import shutil

import pysam
from segment_reads import parse_qname

# functions for initializing new buffer lists
buf_init_funcs = [lambda: (([],),([],)),
                  lambda: (([],[]), ([],[])),
                  lambda: (([],[],[]), ([],[],[])),
                  lambda: (([],[],[],[]), ([],[],[],[])),
                  lambda: (([],[],[],[],[]), ([],[],[],[],[])),
                  lambda: (([],[],[],[],[],[]), ([],[],[],[],[],[])),
                  lambda: (([],[],[],[],[],[],[]), ([],[],[],[],[],[],[])),
                  lambda: (([],[],[],[],[],[],[],[]), ([],[],[],[],[],[],[],[])),
                  lambda: (([],[],[],[],[],[],[],[],[]), ([],[],[],[],[],[],[],[],[])),
                  lambda: (([],[],[],[],[],[],[],[],[],[]), ([],[],[],[],[],[],[],[],[],[]))]
    
def parse_segmented_pe_sam(samfh, maxlen=100000):
    '''
    reads are padded with 3 characters: '_xy', where 
    x=mate number and y=segment number.    
    '''
    qname_ind_map = {}
    start_ind = 0
    end_ind = 0
    cur_ind = 0
    buf = [None] * maxlen
    for read in samfh:
        qname, mate, seg, num_segs = parse_qname(read.qname)
        # set read flags to reflect that this is paired-end
        # data but the reads have not been joined into pairs
        read.qname = qname
        read.is_paired = True
        read.is_proper_pair = False
        read.mate_is_unmapped = True
        read.mrnm = -1
        read.mpos = -1
        read.isize = 0
        if mate == 0:
            read.is_read1 = True
            read.is_read2 = False
        else:
            read.is_read1 = False
            read.is_read2 = True
        # see whether this read is already in the buffer
        if qname not in qname_ind_map:
            buf_size = len(qname_ind_map)
            # test if buffer has become full
            if buf_size == maxlen:
                # buffer full so return first read
                yield buf[start_ind]                
                # delete the read qname to decrease the buffer size by
                # one and allow parser to iterate until it is full again
                return_qname = buf[start_ind][0][0][0].qname
                del qname_ind_map[return_qname]
                # advance start index
                start_ind += 1
                if start_ind == maxlen:
                    start_ind = 0
            # reset end index in buffer            
            qname_ind_map[qname] = end_ind
            buf[end_ind] = buf_init_funcs[num_segs - 1]()
            # get current index for insertion
            cur_ind = end_ind
            # advance end index
            end_ind += 1
            if end_ind == maxlen:
                end_ind = 0
        else:
            # grab buffer index for this read qname
            cur_ind = qname_ind_map[qname]
        # add read to buffer
        buf[cur_ind][mate][seg].append(read)        
    # now empty the rest of the buffer
    buf_size = len(qname_ind_map)
    while buf_size > 0:
        # buffer full so return first read
        yield buf[start_ind]                
        # delete the read qname to decrease the buffer size by
        # one and allow parser to iterate until it is full again
        return_qname = buf[start_ind][0][0][0].qname
        del qname_ind_map[return_qname]
        # advance start index
        start_ind += 1
        if start_ind == maxlen:
            start_ind = 0
        buf_size = len(qname_ind_map)

def build_segment_alignment_dict(aln_dict, reads, seg_num, seg_offset, num_segs):
    for r in reads:
        if r.is_unmapped:
            continue
        # the anchor position is not necessarily the beginning of 
        # the read but is used as a way of checking if adjacent
        # segments are compatible for joining
        if r.is_reverse:
            anchor_pos = r.pos + seg_offset
            seg_ind = num_segs - 1 - seg_num
        else:
            anchor_pos = r.pos - seg_offset
            seg_ind = seg_num
        #anchor_pos = r.pos + seg_offset if r.is_reverse else r.pos - seg_offset            
        assert aln_dict[(r.rname, r.is_reverse, anchor_pos)][seg_ind] == None        
        aln_dict[(r.rname, r.is_reverse, anchor_pos)][seg_ind] = r

def make_unmapped_copy(r):
    a = pysam.AlignedRead()
    a.qname = r.qname
    a.seq = r.seq
    a.qual = r.qual
    a.is_unmapped = True
    a.is_qcfail = False
    a.is_paired = True
    a.is_proper_pair = False
    a.mate_is_unmapped = True
    a.mrnm = -1
    a.mpos = -1
    a.is_read1 = r.is_read1
    a.is_read2 = r.is_read2
    a.isize = 0
    a.mapq = 255
    a.is_reverse = False
    a.rname = -1
    a.pos = 0
    a.cigar = ()
    a.tags = []
    return a

def find_valid_segment_alignments(read_mappings):
    # build a map of reference positions to read segment positions
    num_segs = len(read_mappings)
    aln_dict = collections.defaultdict(lambda: [None] * num_segs)
    unmapped_segs = [None] * num_segs    
    seg_offset = 0
    for seg_num, seg_mappings in enumerate(read_mappings):
        # make an unmapped version of each segment to use in joining
        if len(seg_mappings) == 1 and seg_mappings[0].is_unmapped == True:
            unmapped_segs[seg_num] = seg_mappings[0]
        else:
            unmapped_segs[seg_num] = make_unmapped_copy(seg_mappings[0])
        #print 'SEGMENT', seg_num, 'offset', seg_offset
        build_segment_alignment_dict(aln_dict, seg_mappings, seg_num, seg_offset, num_segs)
        seg_offset += len(seg_mappings[0].seq)

    reads_to_join = []
    if len(aln_dict) == 0:
        # if there are no mappings, then all the segments must be non-mapping
        # and must only have one entry at the 0th position in the segment
        # mapping array.  create an unmapped entry in this case.
        reads_to_join.append([unmapped_segs])
        #seg_reads = [read_mappings[i][0] for i in xrange(num_segs)]
        #reads_to_join.append([seg_reads])
    else:
        # there are some segment mappings, but some of the segments could
        # still be unmapped.  here we create a list of joined sub-segments
        # that together comprise a full read by joining consecutive mapped
        # and unmapped segments
        # find the alignments where the maximum number of segments are
        # joined (the set of best alignments) by decorating the lists with
        # the number of mapped segments then using it as a sort key
        sorted_segment_hits = []
        for mapping_info, segment_hits in aln_dict.iteritems():
            rname, is_reverse, anchor_pos = mapping_info
            num_mapping_segs = 0
            for seg_ind, hit in enumerate(segment_hits):
                # replace missing segments (e.g. 'None') in segment hit list
                # with the unmapped read segment from the main read mappings list
                if hit is None:
                    # find original segment index in read mapping list
                    orig_seg_num = num_segs - 1 - seg_ind if is_reverse else seg_ind
                    # replace with unmapped read
                    segment_hits[seg_ind] = unmapped_segs[orig_seg_num]                    
                    #segment_hits[seg_ind] = read_mappings[orig_seg_num][0]                    
                else:
                    # keep track of number of good mapping segments
                    num_mapping_segs += 1 
            #num_mapping_segs = sum(0 if (i is None) else 1 for i in segment_hits)
            sorted_segment_hits.append((num_mapping_segs, is_reverse, segment_hits)) 
        # get the most segments that mapped contiguously
        sorted_segment_hits.sort(reverse=True)    
        best_num_segs = sorted_segment_hits[0][0]
        #print 'best num segs', best_num_segs
        #print 'hit indexes', [(x[0], x[1], len(x[2])) for x in sorted_segment_hits]
        for num_mapping_segs, is_reverse, segment_hits in sorted_segment_hits:
            if num_mapping_segs < best_num_segs:
                break
            # initialize state for joining
            split_reads = []
            seg_reads = [segment_hits[0]]
            prev_unmapped = seg_reads[-1].is_unmapped
            # find adjacent segments to join
            for seg_num in xrange(1, len(segment_hits)):
                seg_read = segment_hits[seg_num]
                unmapped = seg_read.is_unmapped
                # transition from mapped -> unmapped or unmapped -> mapped
                # requires splitting the read into pieces
                if prev_unmapped != unmapped:
                    split_reads.append(seg_reads)
                    seg_reads = []
                # update unmapped state
                prev_unmapped = unmapped
                seg_reads.append(seg_read)
            # add remaining segments
            if len(seg_reads) > 0:
                split_reads.append(seg_reads)
            # add the split reads to the main list of reads to join
            reads_to_join.append(split_reads)
    return reads_to_join

def join_segments(read_mappings):
    # search for segment matches
    joined_read_hits = find_valid_segment_alignments(read_mappings)
    return joined_read_hits

def parse_MD_tag(val):
    x = 0
    mdops = []
    for y in xrange(len(val)):
        if val[y].isalpha():
            offset = int(val[x:y])
            base = val[y]
            mdops.append(offset)
            mdops.append(base)
            x = y + 1
    if x < len(val):
        mdops.append(int(val[x:]))
    return mdops

def merge_MD_tags(vals):
    mdops = parse_MD_tag(vals[0])
    for val in vals[1:]:
        nextops = parse_MD_tag(val)
        if isinstance(mdops[-1], int) and isinstance(nextops[0], int):
            mdops[-1] += nextops[0]
            nextops = nextops[1:]
        mdops.extend(nextops)
    return ''.join(map(str, mdops))

#def copy_read(r):
#    a = pysam.AlignedRead()
#    a.qname = r.qname
#    a.seq = r.seq
#    a.qual = r.qual
#    a.is_unmapped = r.is_unmapped
#    a.is_qcfail = r.is_qcfail
#    a.is_paired = r.is_paired
#    a.is_proper_pair = r.is_proper_pair
#    a.mate_is_unmapped = r.mate_is_unmapped
#    a.mrnm = r.mrnm
#    a.mpos = r.mpos
#    a.is_read1 = r.is_read1
#    a.is_read2 = r.is_read2
#    a.isize = r.isize
#    a.mapq = r.mapq
#    a.is_reverse = r.is_reverse
#    a.rname = r.rname
#    a.pos = r.pos
#    a.cigar = r.cigar
#    a.tags = r.tags
#    return a
               
def make_joined_read(mate, reads, tags=None):
    if tags is None:
        tags = []
    a = pysam.AlignedRead()
    # create paired-end reads but do not mark them
    # as proper pairs and set all mate information
    # to 'unmapped'
    a.qname = reads[0].qname
    a.seq = ''.join(r.seq for r in reads)
    a.qual = ''.join(r.qual for r in reads)
    a.is_paired = True
    a.is_proper_pair = False
    a.mate_is_unmapped = True
    a.mrnm = -1
    a.mpos = -1
    if mate == 0:
        a.is_read1 = True
        a.is_read2 = False
    else:
        a.is_read1 = False
        a.is_read2 = True
    a.isize = 0
    a.mapq = 255
    a.is_unmapped = reads[0].is_unmapped
    if a.is_unmapped:
        a.rname = -1
        a.pos = 0
    else:
        a.is_reverse = reads[0].is_reverse
        a.rname = reads[0].rname
        a.pos = reads[0].pos
        a.cigar = ((0, len(a.seq)),)
        # compute edit dist
        edit_dist = 0
        for r in reads:
            edit_dist += r.opt('NM')
        tags.append(('NM', edit_dist))
        # compute mismatches to reference (MD)
        tags.append(('MD', merge_MD_tags([r.opt('MD') for r in reads])))
    a.tags = tags
    # TODO: check to see all reads are either mapped or unmapped
    # TODO: check strand assignment
    # TODO: check rname    
    assert all(r.is_unmapped == a.is_unmapped for r in reads)
    assert all(r.is_reverse == a.is_reverse for r in reads)
    assert all(r.rname == a.rname for r in reads)
    return a


def join_segmented_pe_reads(input_sam_file, output_bam_file):
    # setup debugging logging messages
    debug_count = 0
    debug_every = 1e5
    debug_next = debug_every
    # open sam file
    infh = pysam.Samfile(input_sam_file, "r")
    #header = infh.header
    outfh = pysam.Samfile(output_bam_file, "wb", template=infh)
    #outfh = pysam.Samfile("-", "w", template=infh)
    # iterate through paired-end alignments
    logging.info("Processing paired alignments...")
    for segmented_pe_reads in parse_segmented_pe_sam(infh):
        debug_count += 1
        if debug_count == debug_next:
            debug_next += debug_every
            logging.info("Processed %d reads" % debug_count)            
        # get alignments    
        for mate, mate_segs in enumerate(segmented_pe_reads):
            joined_read_mappings = join_segments(mate_segs)            
            # total number of mapping hits
            num_mappings = len(joined_read_mappings)
            #print joined_read_mappings
            for mapping_ind, joined_segments in enumerate(joined_read_mappings):
                # number of split alignment hits
                num_split_hits = len(joined_segments)
                # make SAM record for each segment
                for split_index, seg_reads in enumerate(joined_segments):
                    tags = [('IH', num_mappings),
                            ('HI', mapping_ind),
                            ('XN', num_split_hits),
                            ('XI', split_index)]            
                    r = make_joined_read(mate, seg_reads, tags=tags)
                    outfh.write(r)


if __name__ == '__main__':
    from optparse import OptionParser
    import sys
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = OptionParser("usage: %prog [options] <insam> <outsam>")
    options, args = parser.parse_args()
    input_sam_file = args[0]
    output_bam_file = args[1]
    logging.debug("Joining segmented paired-end mappings")
    logging.debug("Input file: %s" % (input_sam_file))
    logging.debug("Output file: %s" % (output_bam_file))
    join_segmented_pe_reads(input_sam_file, output_bam_file)
