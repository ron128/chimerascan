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
import operator
import collections
import pysam

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
    """
    generator function to parse and return a tuple of 
    lists of reads
    """
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

def group_read_pairs(pe_reads):
    """
    Given tuple of ([read1 reads],[read2 reads]) paired-end read alignments
    return mate-pairs and unpaired reads
    """
    # group paired reads
    paired_reads = ([],[])
    unpaired_reads = ([],[])
    for rnum,reads in enumerate(pe_reads):
        for r in reads:
            if r.is_proper_pair:
                paired_reads[rnum].append(r)
            else:
                unpaired_reads[rnum].append(r)
    # check if we have at least one pair
    pairs = []
    if all((len(reads) > 0) for reads in paired_reads):
        # index read1 by mate reference name and position
        rdict = collections.defaultdict(lambda: collections.deque())
        for r in paired_reads[0]:
            rdict[(r.rnext,r.pnext)].append(r)
        # iterate through read2 and get mate pairs
        for r2 in paired_reads[1]:
            r1 = rdict[(r2.tid,r2.pos)].popleft()
            pairs.append((r1,r2))
    return pairs, unpaired_reads

def select_best_scoring_pairs(pairs):
    """
    return the set of read pairs (provided as a list of tuples) with
    the highest summed alignment score
    """
    if len(pairs) == 0:
        return []
    # gather alignment scores for each pair
    pair_scores = [(pair[0].opt('AS') + pair[1].opt('AS'), pair) for pair in pairs]
    pair_scores.sort(key=operator.itemgetter(0))
    best_score = pair_scores[0][0]
    best_pairs = [pair_scores[0][1]]
    for score,pair in pair_scores[1:]:
        if score < best_score:
            break
        best_pairs.append(pair)
    return best_pairs

def select_primary_alignments(reads):
    """
    return only reads that lack the secondary alignment bit
    """
    if len(reads) == 0:
        return []
    # sort reads by number of mismatches
    unmapped_reads = []
    primary_reads = []
    for r in reads:
        if r.is_unmapped:
            unmapped_reads.append(r)
        elif not r.is_secondary:
            primary_reads.append(r)
    if len(primary_reads) == 0:
        assert len(unmapped_reads) > 0
        return unmapped_reads
    return primary_reads

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
    sorted_reads.extend((worst_nm+1, r) for r in unmapped_reads)
    # choose reads within a certain mismatch tolerance
    best_reads = []
    for mismatches, r in sorted_reads:
        if mismatches > (best_nm + mismatch_tolerance):
            break
        best_reads.append(r)
    return best_reads

def copy_read(r):
    a = pysam.AlignedRead()
    a.qname = r.qname
    a.seq = r.seq
    a.flag = r.flag
    a.tid = r.tid
    a.pos = r.pos
    a.mapq = r.mapq
    a.cigar = r.cigar
    a.rnext = r.rnext
    a.pnext = r.pnext
    a.isize = r.isize
    a.qual = r.qual
    a.tags = list(r.tags)
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
    r1.rnext = r2.tid
    r1.pnext = r2.pos
    tags1 = collections.OrderedDict(r1.tags)
    tags1.update(tags)
    r1.tags = tags1.items()
    # convert read2 to paired-end        
    r2.is_paired = True
    r2.is_proper_pair = True
    r2.is_read2 = True
    r2.mate_is_reverse = r1.is_reverse
    r2.mate_is_unmapped = r1.is_unmapped
    r2.rnext = r1.tid
    r2.pnext = r1.pos
    tags2 = collections.OrderedDict(r2.tags)
    tags2.update(tags)
    r2.tags = tags2.items()
    # compute insert size
    if r1.tid != r2.tid:
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

