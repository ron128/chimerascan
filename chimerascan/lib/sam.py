'''
Created on Jun 2, 2011

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
from chimerascan import pysam
from seq import DNA_reverse_complement

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

def parse_reads_by_qname(samfh):
    """
    generator function to parse and return lists of
    reads that share the same qname
    """    
    reads = []
    for read in samfh:        
        if len(reads) > 0 and read.qname != reads[-1].qname:
            yield reads
            reads = []
        reads.append(read)
    if len(reads) > 0:
        yield reads

def parse_pe_reads(bamfh):
    pe_reads = ([], [])
    # reads must be sorted by qname
    num_reads = 0
    prev_qname = None
    for read in bamfh:
        # get read attributes
        qname = read.qname
        readnum = 1 if read.is_read2 else 0
        # if query name changes we have completely finished
        # the fragment and can reset the read data
        if num_reads > 0 and qname != prev_qname:
            yield pe_reads
            # reset state variables
            pe_reads = ([], [])
            num_reads = 0
        pe_reads[readnum].append(read)
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

def copy_read(r):
    a = pysam.AlignedRead()
    a.qname = r.qname
    a.seq = r.seq
    a.flag = r.flag
    a.rname = r.rname
    a.pos = r.pos
    a.mapq = r.mapq
    a.cigar = r.cigar
    a.mrnm = r.mrnm
    a.mpos = r.mpos
    a.isize = r.isize
    a.qual = r.qual
    a.tags = r.tags
    return a

def soft_pad_read(fq, r):
    """
    'fq' is the fastq record
    'r' in the AlignedRead SAM read
    """    
    # make sequence soft clipped
    ext_length = len(fq.seq) - len(r.seq)
    cigar_softclip = [(CIGAR_S, ext_length)]
    cigar = r.cigar
    # reconstitute full length sequence in read
    if r.is_reverse:
        seq = DNA_reverse_complement(fq.seq)
        qual = fq.qual[::-1]        
        if (cigar is not None) and (ext_length > 0):
            cigar = cigar_softclip + cigar
    else:
        seq = fq.seq
        qual = fq.qual
        if (cigar is not None) and (ext_length > 0):
            cigar = cigar + cigar_softclip
    # replace read field
    r.seq = seq
    r.qual = qual
    r.cigar = cigar

def pair_reads(r1, r2, tags=None):
    '''
    fill in paired-end fields in SAM record
    '''
    if tags is None:
        tags = []
    # convert read1 to paired-end
    r1.is_paired = True
    r1.is_proper_pair = True
    r1.is_read1 = True
    r1.mate_is_reverse = r2.is_reverse
    r1.mate_is_unmapped = r2.is_unmapped
    r1.mpos = r2.pos
    r1.mrnm = r2.rname
    r1.tags = r1.tags + tags
    # convert read2 to paired-end        
    r2.is_paired = True
    r2.is_proper_pair = True
    r2.is_read2 = True
    r2.mate_is_reverse = r1.is_reverse
    r2.mate_is_unmapped = r1.is_unmapped
    r2.mpos = r1.pos
    r2.mrnm = r1.rname
    r2.tags = r2.tags + tags
    # compute insert size
    if r1.rname != r2.rname:
        r1.isize = 0
        r2.isize = 0
    elif r1.pos > r2.pos:
        isize = r1.aend - r2.pos
        r1.isize = -isize
        r2.isize = isize
    else:
        isize = r2.aend - r1.pos
        r1.isize = isize
        r2.isize = -isize

def get_clipped_interval(r):
    cigar = r.cigar
    padstart, padend = r.pos, r.aend
    if len(cigar) > 1:
        if (cigar[0][0] == CIGAR_S or
            cigar[0][0] == CIGAR_H):
            padstart -= cigar[0][1]
        elif (cigar[-1][0] == CIGAR_S or
            cigar[-1][0] == CIGAR_H):
            padend += cigar[-1][1]
    return padstart, padend

