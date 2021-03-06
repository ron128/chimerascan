'''
Created on Jan 23, 2011

@author: mkiyer

chimerascan: chimeric transcript discovery using RNA-seq

Copyright (C) 2011 Matthew Iyer

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''
import re
import collections

from chimerascan.lib.seq import parse_segmented_qname

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

def parse_fastq_qname_mate(line_iter):
    mate_re = re.compile(r'/(\d)$')    
    try:
        while True:
            # qname
            qname = line_iter.next().rstrip()[1:]
            qname,mate = mate_re.split(qname)[0:2]
            # skip 3 lines
            line_iter.next()
            line_iter.next()
            line_iter.next()
            yield qname, int(mate) - 1
    except StopIteration:
        pass

def fix_pe_alignment_ordering(samfh, fastq_iter, is_paired=True, maxlen=100000):
    num_mates = 2 if is_paired else 1    
    # function for initializing new buffer list
    buf_init_func = lambda: tuple(list() for m in xrange(num_mates))
    # initialize the qname dictionary to match the fastq file    
    buf = collections.deque()
    qname_read_dict = {}
    qname_iter = parse_fastq_qname(fastq_iter)
    for read in samfh:
        if is_paired:
            mate = 0 if read.is_read1 else 1
        else:
            mate = 0
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

def fix_sr_alignment_ordering(samfh, fastq_iter, maxlen=100000):
    # initialize the qname dictionary to match the fastq file    
    buf = collections.deque()
    qname_read_dict = {}
    qname_iter = parse_fastq_qname_mate(fastq_iter)    
    qname_mate_re = re.compile(r'/(\d)$')    
    for read in samfh:
        # get mate from SAM file
        read_qname, read_mate = qname_mate_re.split(read.qname)[0:2]
        read_mate = int(read_mate) - 1
        # set key for indexing reads
        key = (read_qname, read_mate)
        # set flags
        read.qname = read_qname
        if read_mate == 0:
            read.is_read1 = True
        elif read_mate == 1:
            read.is_read2 = True
        else:
            assert False
        # check if this read is already in the buffer
        if key not in qname_read_dict:
            # if buffer full empty the first entries
            while len(buf) >= maxlen:
                # get first qname in buf
                first_key = buf.popleft()
                # return reads at this qname, then delete them
                yield qname_read_dict[first_key]
                del qname_read_dict[first_key]
            # add new qnames to buffer
            while True:                
                # get next qname from fastq file and add it to the queue
                next_key = qname_iter.next()
                buf.append(next_key)
                qname_read_dict[next_key] = list()
                # if the next qname in the fastq file is the same as the
                # read qname, then we can exit the loop
                if next_key == key:
                    break
        # add read to buffer
        qname_read_dict[key].append(read)
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
        qname, mate, seg, num_segs = parse_segmented_qname(read.qname)
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
        if qname not in qname_read_dict:
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
                if next_qname not in qname_read_dict:
                    # initialize qname entry
                    buf.append(next_qname)
                    qname_read_dict[next_qname] = buf_init_func(num_segs)
                # if the next qname in the fastq file is the same as the
                # read qname, then we can exit the loop
                if next_qname == read.qname:
                    break
        #print 'READ1', read.is_read1, 'QNAME BEFORE', read.qname, qname_read_dict[read.qname]  
        # add read to buffer
        qname_read_dict[qname][mate][seg].append(read)
        #print 'QNAME AFTER', qname, qname_read_dict[read.qname]
    # empty remaining entries in buffer
    while len(buf) > 0:
        yield qname_read_dict[buf.popleft()]
