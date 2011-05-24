'''
Created on Apr 29, 2011

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
import array
import logging

from base import NO_STRAND

def parse_pe_reads(bamfh):
    pe_reads = ([], [])
    # reads must be binned by qname, mate, hit, and segment
    # so initialize to mate 0, hit 0, segment 0
    num_reads = 0
    prev_qname = None
    for read in bamfh:
        # get read attributes
        qname = read.qname
        if read.is_read1:
            mate = 0
        elif read.is_read2:
            mate = 1
        # if query name changes we have completely finished
        # the fragment and can reset the read data
        if num_reads > 0 and qname != prev_qname:
            yield pe_reads
            # reset state variables
            pe_reads = ([], [])
            num_reads = 0
        pe_reads[mate].append(read)
        prev_qname = qname
        num_reads += 1
    if num_reads > 0:
        yield pe_reads

CIGAR_M = 0 #match  Alignment match (can be a sequence match or mismatch)
CIGAR_I = 1 #insertion  Insertion to the reference
CIGAR_D = 2 #deletion  Deletion from the reference
CIGAR_N = 3 #skip  Skipped region from the reference
CIGAR_S = 4 #softclip  Soft clip on the read (clipped sequence present in <seq>)
CIGAR_H = 5 #hardclip  Hard clip on the read (clipped sequence NOT present in <seq>)
CIGAR_P = 6 #padding  Padding (silent deletion from the padded reference sequence)

def get_genomic_intervals(read):
    intervals = []
    rseq = read.seq
    qseq = array.array('c')
    qstart = 0
    astart = read.pos
    aend = astart
    for op,length in read.cigar:
        if (op == CIGAR_D):
            aend += length
        elif (op == CIGAR_I) or (op == CIGAR_S):
            qstart += length
        elif (op == CIGAR_M):            
            qseq.fromstring(rseq[qstart:qstart + length])
            qstart += length
            aend += length
        elif (op == CIGAR_N):
            if aend > astart:
                if len(qseq) != (aend - astart):
                    logging.error("Read %s has aend != astart" % (str(read)))
                else:
                    intervals.append((astart, aend, qseq))
            astart = aend + length
            aend = astart
            qseq = array.array('c')
    if aend > astart:
        if len(qseq) != (aend - astart):
            logging.error("Read %s has aend != astart" % (str(read)))
        else:
            intervals.append((astart, aend, qseq))
    if aend != read.aend:
        logging.error("Read %s has aend != read.aend" % (str(read)))
    return intervals


def get_insert_size(read1, read2):
    # compute the total span of the reads
    if read2.pos < read1.pos:
        span = read1.aend - read2.pos
    else:
        span = read2.aend - read1.pos
    # subtract any "gaps" or "skips" in the alignment
    skips = 0
    for r in (read1, read2):
        for op,length in r.cigar:
            if (op == CIGAR_D):
                skips += length
            elif (op == CIGAR_I) or (op == CIGAR_S):
                span += length
            elif (op == CIGAR_N):
                skips += length
    return span - skips

def get_strand(read, tag="XS"):
    try:
        return read.opt(tag)
    except KeyError:
        return NO_STRAND
