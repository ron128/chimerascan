'''
Created on Jan 23, 2011

@author: mkiyer
'''
import re
import collections

from segment_reads import parse_qname

def parse_fastq_qname(line_iter):
    try:        
        qname = line_iter.next().rstrip()[1:]
        newqname = re.split(r'/\d$', qname)[0]
        suffix_length = len(qname) - len(newqname)                    
        # skip 3 lines
        line_iter.next()
        line_iter.next()
        line_iter.next()
        yield newqname
        while True:
            # qname
            qname = line_iter.next().rstrip()[1:]
            qname = qname[:len(qname)-suffix_length]
            # skip 3 lines
            line_iter.next()
            line_iter.next()
            line_iter.next()
            yield qname
    except StopIteration:
        pass

def fix_alignment_ordering(samfh, fastq_iter, is_paired=True, maxlen=100000):
    num_mates = 2 if is_paired else 1    
    # function for initializing new buffer list
    buf_init_func = lambda: tuple(list() for m in xrange(num_mates))
    # initialize the qname dictionary to match the fastq file    
    buf = collections.deque()
    qname_read_dict = {}
    qname_iter = parse_fastq_qname(fastq_iter)
    for read in samfh:
        mate = 0 if read.is_read1 else 1
        # check if this read is already in the buffer
        if read.qname not in qname_read_dict:
            # if buffer full empty the first entries
            while len(buf) >= maxlen:
                # get first qname in buf
                first_qname = buf.popleft()
                # return reads at this qname, then delete them
                yield qname_read_dict[first_qname]
                del qname_read_dict[first_qname]
            # add new qnames to buffer
            while True:                
                # get next qname from fastq file and add it to the queue
                next_qname = qname_iter.next()
                buf.append(next_qname)
                qname_read_dict[next_qname] = buf_init_func()                
                # if the next qname in the fastq file is the same as the
                # read qname, then we can exit the loop
                if next_qname == read.qname:
                    break
        # add read to buffer
        qname_read_dict[read.qname][mate].append(read)
    # empty remaining entries in buffer
    while len(buf) > 0:
        yield qname_read_dict[buf.popleft()]


def fix_segmented_alignment_ordering(samfh, fastq_iter, is_paired=True, maxlen=100000):
    num_mates = 2 if is_paired else 1    
    # function for initializing new buffer list
    buf_init_func = lambda num_segs: tuple(tuple(list() for x in xrange(num_segs)) 
                                           for m in xrange(num_mates))
    # initialize the qname dictionary to match the fastq file
    buf = collections.deque()
    qname_read_dict = {}
    qname_iter = parse_fastq_qname(fastq_iter)
    for read in samfh:
        qname, mate, seg, num_segs = parse_qname(read.qname)
        # set read flags to reflect that this is paired-end
        # data but the reads have not been joined into pairs
        read.qname = qname
        read.is_paired = is_paired
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
        # check if this read is already in the buffer
        if read.qname not in qname_read_dict:
            # if buffer full empty the first entries
            while len(buf) >= maxlen:
                # get first qname in buf
                first_qname = buf.popleft()
                # return reads at this qname, then delete them
                yield qname_read_dict[first_qname]
                del qname_read_dict[first_qname]
            # add new qnames to buffer
            while True:
                # get next qname from fastq file and add it to the queue
                next_qname = qname_iter.next()
                buf.append(next_qname)
                # initialize qname entry
                qname_read_dict[next_qname] = buf_init_func(num_segs)
                # if the next qname in the fastq file is the same as the
                # read qname, then we can exit the loop
                if next_qname == read.qname:
                    break
        # add read to buffer
        qname_read_dict[read.qname][mate][seg].append(read)
    # empty remaining entries in buffer
    while len(buf) > 0:
        yield qname_read_dict[buf.popleft()]



