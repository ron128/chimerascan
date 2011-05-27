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
import operator

from base import NO_STRAND

#
# constants used for CIGAR alignments
#
CIGAR_M = 0 #match  Alignment match (can be a sequence match or mismatch)
CIGAR_I = 1 #insertion  Insertion to the reference
CIGAR_D = 2 #deletion  Deletion from the reference
CIGAR_N = 3 #skip  Skipped region from the reference
CIGAR_S = 4 #softclip  Soft clip on the read (clipped sequence present in <seq>)
CIGAR_H = 5 #hardclip  Hard clip on the read (clipped sequence NOT present in <seq>)
CIGAR_P = 6 #padding  Padding (silent deletion from the padded reference sequence)

# custom read tags
class SamTags:
    RTAG_NUM_PARTITIONS = "XP"
    RTAG_PARTITION_IND = "XH"
    RTAG_NUM_SPLITS = "XN"
    RTAG_SPLIT_IND = "XI"
    RTAG_NUM_MAPPINGS = "IH"
    RTAG_MAPPING_IND = "HI"
    RTAG_BOWTIE_MULTIMAP = "XM"


def parse_pe_reads(bamfh):
    pe_reads = ([], [])
    # reads must be sorted by qname
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

def parse_unpaired_pe_reads(bamfh):
    """
    parses alignments that were aligned in single read mode
    and hence all hits are labeled as 'read1' and lack mate
    information.  instead the read1 read2 information is
    attached to the 'qname' field
    """
    pe_reads = ([], [])
    num_reads = 0
    prev_qname = None
    for read in bamfh:
        # extract read1/2 from qname
        readnum = int(read.qname[-1])
        if readnum == 1:
            read.is_read1 = True
            mate = 0
        elif readnum == 2:
            mate = 1
            read.is_read2 = True
        # reconstitute correct qname
        qname = read.qname[:-2]
        read.qname = qname
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

def select_best_mismatch_strata(reads, mismatch_tolerance=0):
    if len(reads) == 0:
        return []
    # sort reads by number of mismatches
    mapped_reads = []
    unmapped_reads = []
    for r in reads:
        if r.is_unmapped:
            unmapped_reads.append(r)
        else:
            mapped_reads.append((r.opt('NM'), r))
    if len(mapped_reads) == 0:
        return unmapped_reads
    sorted_reads = sorted(mapped_reads, key=operator.itemgetter(0))
    best_nm = sorted_reads[0][0]
    worst_nm = sorted_reads[-1][0]
    sorted_reads.extend((worst_nm, r) for r in unmapped_reads)
    # choose reads within a certain mismatch tolerance
    best_reads = []
    for mismatches, r in sorted_reads:
        if mismatches > best_nm + mismatch_tolerance:
            break
        best_reads.append(r)
    return best_reads

def parse_multihit_alignments(samfh):
    buf = []
    ind = 0
    for read in samfh:
        if (ind > 0) and (read.qname != buf[ind-1].qname):
            yield buf[:ind]
            ind = 0
        if ind < len(buf):
            buf[ind] = read
        else:
            buf.append(read)
        ind += 1
    if ind > 0:
        yield buf[:ind]

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