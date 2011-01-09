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

def write_discordant_reads(reads, num_hits, outfh):
    multihit = False
    for r in reads:
        # do not want more than one unmapped read per mate
        if r.is_unmapped and multihit:
            continue
        multihit = True
        #logging.debug("read: %s num_hits: %d" % (r.qname, num_hits))
        r.tags = r.tags + [('NH', num_hits)]
        outfh.write(r)

CIGAR_M = 0 #match  Alignment match (can be a sequence match or mismatch)
CIGAR_I = 1 #insertion  Insertion to the reference
CIGAR_D = 2 #deletion  Deletion from the reference
CIGAR_N = 3 #skip  Skipped region from the reference
CIGAR_S = 4 #softclip  Soft clip on the read (clipped sequence present in <seq>)
CIGAR_H = 5 #hardclip  Hard clip on the read (clipped sequence NOT present in <seq>)
CIGAR_P = 6 #padding  Padding (silent deletion from the padded reference sequence)

def get_aligned_read_intervals(read):
    intervals = []
    # insert read into cluster tree
    astart,aend = read.pos, read.pos
    for op,length in read.cigar:
        if length == 0: continue
        if (op == CIGAR_I) or (op == CIGAR_S) or (op == CIGAR_H): continue
        if (op == CIGAR_P): assert False 
        if (op == CIGAR_N):
            assert astart != aend
            intervals.append((astart, aend))
            #print read.qname, read.cigar, ref, astart, aend
            astart = aend + length
        aend += length
    assert astart != aend
    if aend > astart:
        #print read.qname, read.cigar, ref, astart, aend
        intervals.append((astart, aend))
    assert aend == read.aend
    return intervals

def partition_mappings_by_strand(read_mappings):
    # bin segments by strand    
    num_segs = len(read_mappings)
    strand_segs = ([list() for x in xrange(num_segs)],
                   [list() for x in xrange(num_segs)])                   
    for seg_num, seg_mappings in enumerate(read_mappings):
        for r in seg_mappings:
            if r.is_unmapped:
                # add unmapped reads to the lists on both
                # strands since the strands are joined
                # separately
                strand_segs[0][seg_num].append(r)
                strand_segs[1][num_segs - 1 - seg_num].append(r)
            else:
                # flip negative strand segments
                if r.is_reverse:
                    seg_ind = num_segs - 1 - seg_num
                else:
                    seg_ind = seg_num
                strand_segs[r.is_reverse][seg_ind].append(r)
    return strand_segs

def build_segment_alignment_dict(aln_dict, reads, seg_num, seg_offset, num_segs):
    for i,r in enumerate(reads):
        if r.is_unmapped:
            continue
        # the anchor position is not necessarily the beginning of 
        # the read but is used as a way of checking if adjacent
        # segments are compatible for joining
        anchor_pos = r.pos + seg_offset if r.is_reverse else r.pos - seg_offset            
        assert aln_dict[(r.rname, r.is_reverse, anchor_pos)][seg_num] == None        
        aln_dict[(r.rname, r.is_reverse, anchor_pos)][seg_num] = i

def find_valid_segment_alignments(read_mappings):
    # build a map of reference positions to read segment positions
    num_segs = len(read_mappings)
    aln_dict = collections.defaultdict(lambda: [None] * num_segs)
    seg_offset = 0
    for seg_num, seg_mappings in enumerate(read_mappings):
        #print 'SEGMENT', seg_num, 'offset', seg_offset
        build_segment_alignment_dict(aln_dict, seg_mappings, seg_num, seg_offset, num_segs)
        seg_offset += len(seg_mappings[0].seq)

    # find the alignments where the maximum number of segments are
    # joined (the set of best alignments) by decorating the lists with
    # the number of mapped segments then using it as a sort key
    decorated_hit_indexes = []
    #for mapping_info, segment_hit_indexes in aln_dict.iteritems():        
    for mapping_info, segment_hit_indexes in aln_dict.iteritems():
        rname, is_reverse, anchor_pos = mapping_info        
        num_mapping_segs = sum(0 if (i is None) else 1 for i in segment_hit_indexes)
        decorated_hit_indexes.append((num_mapping_segs, is_reverse, segment_hit_indexes)) 
    decorated_hit_indexes.sort(reverse=True)    

    # search reference map and join reads where possible
    reads_to_join = []
    if len(decorated_hit_indexes) == 0:
        # if there are no mappings, then all the segments must be non-mapping
        # and must only have one entry at the 0th position in the segment
        # mapping array.  create an unmapped entry in this case.
        seg_reads = [read_mappings[i][0] for i in xrange(num_segs)]
        reads_to_join.append([seg_reads])
    else:
        # there are some segment mappings, but some of the segments could
        # still be unmapped.  here we create a list of joined sub-segments
        # that together comprise a full read by joining consecutive mapped
        # and unmapped segments
        best_num_segs = decorated_hit_indexes[0][0]
        #print 'best num segs', best_num_segs
        #print 'hit indexes', decorated_hit_indexes
        for num_mapping_segs, is_reverse, segment_hit_indexes in decorated_hit_indexes:
            if num_mapping_segs < best_num_segs:
                break
            # initialize state for joining
            split_reads = []
            ind = segment_hit_indexes[0]
            prev_unmapped = (ind == None)
            correct_ind = 0 if prev_unmapped else ind
            seg_reads = [read_mappings[0][correct_ind]]        
            # find adjacent segments to join
            for seg_num in xrange(1, len(segment_hit_indexes)):
                ind = segment_hit_indexes[seg_num]
                unmapped = (ind == None)
                correct_ind = 0 if unmapped else ind
                # transition from mapped -> unmapped or unmapped -> mapped
                # requires splitting the read into pieces
                if prev_unmapped != unmapped:
                    split_reads.append(seg_reads)
                    seg_reads = []
                # update unmapped state
                prev_unmapped = unmapped
                seg_reads.append(read_mappings[seg_num][correct_ind])                
            # add remaining segments
            if len(seg_reads) > 0:
                split_reads.append(seg_reads)
            # add the split reads to the main list of reads to join
            # TODO: is this the best way to join strands?
            if is_reverse:
                split_reads.reverse()
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
            mdops.append((offset, base))
            x = y + 1
#    if x < len(val):
#        mdops.append(int(val[x:]))
    return mdops

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
    a.tags = tags
    #a.tags = ( ("NM", 1),
    #           ("RG", "L1") )
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
    #outfh = pysam.Samfile(output_bam_file, "wb", template=infh)
    outfh = pysam.Samfile("-", "w", template=infh)
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
