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

ReorderBufferItem = collections.namedtuple('ReorderBufferItem', ("fqrec", "reads"))

def fix_alignment_ordering(samfh, fqiters, 
                           pe_sr_mode=False,
                           maxlen=100000):
    # function for initializing new buffer entry
    buf_init_func = lambda fqrecs: tuple(ReorderBufferItem(fq, []) for fq in fqrecs)
    # initialize the qname dictionary to match the fastq file    
    buf = collections.deque()
    qname_read_dict = {}
    qname_mate_re = re.compile(r'/(\d)$')
    for read in samfh:
        # PE-SR mode means that the reads were paired in sequencing
        # but aligned separately.  The function uses the /1 and /2
        # suffixes in the reads to join them during buffer reordering
        if pe_sr_mode:
            # get read num (1 or 2) from the qname field of SAM read
            read_qname, readnum = qname_mate_re.split(read.qname)[0:2]
            readnum = int(readnum) - 1
            # set flags
            read.is_paired = True
            read.qname = read_qname
            if readnum == 0:
                read.is_read1 = True
            elif readnum == 1:
                read.is_read2 = True
            else:
                assert False
        # if not PE-SR mode then we can trust the 'is_read1' and 'is_read2'
        # attributes of the SAM read
        else:
            if read.is_read2:
                readnum = 1
            else:
                readnum = 0
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
                fqrecs = [it.next() for it in fqiters]
                next_qname = fqrecs[0].qname
                buf.append(next_qname)
                qname_read_dict[next_qname] = buf_init_func(fqrecs)
                # if the next qname in the fastq file is the same as the
                # read qname, then we can exit the loop
                if next_qname == read.qname:
                    break
        # add read to buffer
        qname_read_dict[read.qname][readnum].reads.append(read)
    # empty remaining entries in buffer
    while len(buf) > 0:
        yield qname_read_dict[buf.popleft()]

# TODO: remove this
#
#def fix_pe_sr_alignment_ordering(samfh, fastq_line_iters, maxlen=100000):
#    # function for initializing new buffer entry
#    buf_init_func = lambda fqrecs: tuple(ReorderBufferItem(fq, []) for fq in fqrecs)
#    # initialize the qname dictionary to match the fastq file    
#    buf = collections.deque()
#    qname_read_dict = {}
#    qname_mate_re = re.compile(r'/(\d)$')
#    fastq_rec_iters = [parse_fastq_record(it) for it in fastq_line_iters]
#    for read in samfh:
#        # get mate from SAM file
#        read_qname, readnum = qname_mate_re.split(read.qname)[0:2]
#        readnum = int(readnum) - 1
#        # set flags
#        read.qname = read_qname
#        if readnum == 0:
#            read.is_read1 = True
#        elif readnum == 1:
#            read.is_read2 = True
#        else:
#            assert False
#        # check if this read is already in the buffer
#        if read.qname not in qname_read_dict:
#            # if buffer full empty the first entries
#            while len(buf) >= maxlen:
#                # get first qname in buf
#                first_qname = buf.popleft()
#                # return reads at this qname, then delete them
#                yield qname_read_dict[first_qname]
#                del qname_read_dict[first_qname]
#            # add new qnames to buffer
#            while True:                
#                # get next qname from fastq file and add it to the queue
#                fqrecs = [it.next() for it in fastq_rec_iters]
#                next_qname = fqrecs[0].qname
#                buf.append(next_qname)
#                qname_read_dict[next_qname] = buf_init_func(fqrecs)
#                # if the next qname in the fastq file is the same as the
#                # read qname, then we can exit the loop
#                if next_qname == read.qname:
#                    break
#        # add read to buffer
#        qname_read_dict[read.qname][readnum].append(read)
#    # empty remaining entries in buffer
#    while len(buf) > 0:
#        yield qname_read_dict[buf.popleft()]
#
#
## TODO: remove or modify
#def fix_sr_alignment_ordering(samfh, fastq_iter, maxlen=100000):
#    pass
## TODO: remove or modify
##    # initialize the qname dictionary to match the fastq file    
##    buf = collections.deque()
##    qname_read_dict = {}
##    qname_iter = parse_fastq_qname_mate(fastq_iter)    
##    qname_mate_re = re.compile(r'/(\d)$')
##
##    for read in samfh:
##        # get mate from SAM file
##        read_qname, read_mate = qname_mate_re.split(read.qname)[0:2]
##        read_mate = int(read_mate) - 1
##        # set key for indexing reads
##        key = (read_qname, read_mate)
##        # set flags
##        read.qname = read_qname
##        if read_mate == 0:
##            read.is_read1 = True
##        elif read_mate == 1:
##            read.is_read2 = True
##        else:
##            assert False
##        # check if this read is already in the buffer
##        if key not in qname_read_dict:
##            # if buffer full empty the first entries
##            while len(buf) >= maxlen:
##                # get first qname in buf
##                first_key = buf.popleft()
##                # return reads at this qname, then delete them
##                yield qname_read_dict[first_key]
##                del qname_read_dict[first_key]
##            # add new qnames to buffer
##            while True:                
##                # get next qname from fastq file and add it to the queue
##                next_key = qname_iter.next()
##                buf.append(next_key)
##                qname_read_dict[next_key] = list()
##                # if the next qname in the fastq file is the same as the
##                # read qname, then we can exit the loop
##                if next_key == key:
##                    break
##        # add read to buffer
##        qname_read_dict[key].append(read)
##    # empty remaining entries in buffer
##    while len(buf) > 0:
##        yield qname_read_dict[buf.popleft()]

