'''
Created on Jan 23, 2011

@author: mkiyer
'''
import re
import collections

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
    buf = collections.deque() 
    # initialize the qname dictionary to match the fastq file
    qname_iter = parse_fastq_qname(fastq_iter)
    qname_ind_map = {}
    for read in samfh:
        mate = 0 if read.is_read1 else 1
        # check if this read is already in the buffer
        if read.qname not in qname_ind_map:
            # if buffer full empty the first entries
            while len(buf) >= maxlen:
                # return first read in buffer
                del qname_ind_map[buf[0][0][0].qname]
                yield buf.popleft()
            # add new qnames to buffer
            while True:
                # get next qname from fastq file and add it to the end of 
                # the buffer                
                next_qname = qname_iter.next()
                next_ind = len(buf)
                qname_ind_map[next_qname] = next_ind
                buf.append(buf_init_func())
                # if the next qname in the fastq file is the same as the
                # read qname, then we can exit the loop
                if next_qname == read.qname:
                    # get current index for insertion
                    cur_ind = next_ind
                    break
        else:
            # grab buffer index for this read qname
            cur_ind = qname_ind_map[read.qname]
        # add read to buffer
        buf[cur_ind][mate].append(read)        
    # empty remaining entries in buffer
    while len(buf) >= maxlen:
        # return first read in buffer
        #del qname_ind_map[buf[0][0][0].qname]
        yield buf.popleft()
