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
import collections
import array
import logging
import random

from chimerascan.bx.intersection import Interval, IntervalTree

# local imports
from sam import parse_pe_reads, CIGAR_N, CIGAR_S, CIGAR_H, CIGAR_P

# SAM CIGAR flags that indicate skipping, padding, or clipping
SKIP_CIGAR_FLAGS = set((CIGAR_N, CIGAR_S, CIGAR_H, CIGAR_P)) 

def build_exon_trees(genes):
    trees = collections.defaultdict(lambda: IntervalTree())
    for g in genes:        
        for e in g.exons:
            start, end = e
            trees[g.chrom].insert_interval(Interval(start, end, strand=g.strand))
    return trees

def find_unambiguous_exon_intervals(genes):
    """
    returns (chrom, start, end, strand) tuples for exon
    intervals that are unique and have no overlapping
    transcripts or exons.    
    """
    trees = build_exon_trees(genes)    
    for g in genes:
        for start,end in g.exons:
            hits = set((hit.start, hit.end, hit.strand) 
                       for hit in trees[g.chrom].find(start, end))
            hits.add((start, end, g.strand))
            if len(hits) == 1:
                yield g.chrom, start, end, g.strand

def sample_fragment_sizes(bamfh, genes, min_isize, max_isize):
    """
    sample fragment size distribution at genes with exons
    larger than the maximum insert size
    """
    # find all exons that are larger than the maximum estimated fragment size
    exons = set(coord for coord in find_unambiguous_exon_intervals(genes)
                if (coord[2] - coord[1]) >= max_isize)
    logging.info("Found %d exons larger than %d" % (len(exons), max_isize))
    refs = set(bamfh.references)
    # stats
    num_reads = 0
    unmapped = 0
    ambiguous = 0
    spliced = 0
    outside_range = 0
    count = 0
    # fetch reads from BAM file at large exons
    for chrom,start,end,strand in exons:
        if chrom not in refs:
            logging.warning("Skipping exon from reference %s not in BAM" % (chrom))
            continue 
        qname_dict = collections.defaultdict(lambda: [])
        for r in bamfh.fetch(chrom, start, end):
            num_reads += 1
            # ignore unmapped reads, qc fail reads, or unpaired reads
            if r.is_unmapped or r.is_qcfail or (not r.is_proper_pair):
                unmapped += 1
                continue
            # ignore multi-mapping reads
            if r.opt('NH') > 1:
                ambiguous += 1
                continue
            # ignore spliced reads
            has_skip = any(x[0] in SKIP_CIGAR_FLAGS for x in r.cigar)
            if has_skip:
                spliced += 1
                continue            
            # group paired-end reads by read name
            qname_dict[r.qname].append(abs(r.isize))
        # keep paired reads with both mates in region
        for isizes in qname_dict.itervalues():
            isizes = set(abs(x) for x in isizes)
            assert len(isizes) == 1
            isize = isizes.pop()
            if (min_isize <= isize <= max_isize):
                count += 1
                yield isize
            else:
                outside_range += 1
    logging.debug("Processed reads: %d" % (num_reads))
    logging.debug("Unique paired frags: %d" % (count))
    logging.debug("Outside range: %d" % (outside_range))
    logging.debug("Unmapped or invalid: %d" % (unmapped))
    logging.debug("Ambiguous: %d" % (ambiguous))
    logging.debug("Spliced or padded: %d" % (spliced))

class InsertSizeDistribution(object):
    
    def __init__(self):
        self.min_isize = None
        self.max_isize = None
        self.arr = None

    def isize_at_percentile(self, per):
        n = sum(self.arr)
        per_n = n * per / 100.0
        count = 0
        for isize,x in enumerate(self.arr): 
            count += x
            if (count >= per_n):
                break
        return isize + self.min_isize
    
    def percentile_at_isize(self, isize):
        if isize < self.min_isize:
            return 0.0
        elif isize > self.max_isize:
            return 100.0
        ind = isize - self.min_isize
        count_le = sum(self.arr[:ind+1])        
        per = 100.0 * count_le / float(sum(self.arr))
        return per

    @property
    def n(self):
        if self.arr is None: return 0
        return sum(self.arr)
    
    def mode(self):
        return self.arr.index(max(self.arr)) + self.min_isize

    def mean(self):
        count = 0
        n = 0        
        for i,x in enumerate(self.arr): 
            count += i*x
            n += x
        if n == 0:
            return None            
        return self.min_isize + (count / float(n))
    
    def std(self):
        # first get mean (not shifted)
        mean = self.mean() - self.min_isize
        if mean is None:
            return None
        n = 0
        std = 0
        for i,x in enumerate(self.arr):
            std = std + x*((i - mean)**2)
            n += x
        std = (std / float(n-1))**0.5
        return std

    def to_file(self, fileh):
        print >>fileh, '\t'.join(["#insert_size", "num_samples"])
        for i,x in enumerate(self.arr):
            print >>fileh, '\t'.join([str(i + self.min_isize), str(x)])        

    @staticmethod
    def from_file(fileh):
        isizes = []
        counts = []
        for line in fileh:
            if line.startswith("#"):
                continue
            fields = line.strip().split('\t')
            i,x = map(int, fields[0:2])
            isizes.append(i)
            counts.append(x)
        d = InsertSizeDistribution()
        d.min_isize = isizes[0]
        d.max_isize = isizes[-1]
        d.arr = array.array('L', counts) 
        return d

    @staticmethod
    def from_random(mean, stdev, min_isize, max_isize, samples=100000):
        """
        initialize from a random sample using normal distribution with 
        mean 'mean' and stdev 'stdev'
        """
        d = InsertSizeDistribution()
        # implement simple checks
        assert min_isize < mean < max_isize
        assert stdev < (max_isize - min_isize)
        # initialize
        d.min_isize = min_isize
        d.max_isize = max_isize
        d.arr = array.array('L', (0 for x in xrange(min_isize, max_isize+1)))
        count = 0
        outside_range = 0
        while True:
            if count > samples:
                break
            isize = int(round(random.normalvariate(mean, stdev),0))
            if (min_isize <= isize <= max_isize):
                # store in array
                d.arr[isize - min_isize] += 1
                count += 1
            else:
                outside_range += 1
        return d

    @staticmethod
    def from_bam(bamfh, min_isize, max_isize, max_samples=None):
        # initialize
        d = InsertSizeDistribution()
        d.min_isize = min_isize
        d.max_isize = max_isize
        d.arr = array.array('L', (0 for x in xrange(min_isize, max_isize+1)))     
        frags = 0   
        count = 0
        outside_range = 0
        unmapped = 0
        isoforms = 0
        for pe_reads in parse_pe_reads(bamfh):
            frags += 1
            if (max_samples is not None) and (count > max_samples):
                break
            # only allow mappings where there is a single
            # insert size (multiple isoforms are ambiguous)
            isizes = set()        
            for r in pe_reads[0]:
                if r.is_unmapped:
                    continue
                # get insert size
                isize = r.isize
                if isize < 0: isize = -isize
                isizes.add(isize)
            # insert size must be within range
            if len(isizes) == 0:
                unmapped += 1
            elif len(isizes) > 1:
                isoforms += 1
            else:
                isize = isizes.pop()
                if (min_isize <= isize <= max_isize):
                    # store in array
                    d.arr[isize - min_isize] += 1
                    count += 1
                else:
                    outside_range += 1
        logging.debug("Processed fragments: %d" % (frags))
        logging.debug("Unique paired frags: %d" % (count))
        logging.debug("Unmapped: %d" % (unmapped))
        logging.debug("Ambiguous (isoforms): %d" % (isoforms))
        logging.debug("Outside range: %d" % (outside_range))
        return d
    
    @staticmethod
    def from_genome_bam(bamfh, genes, min_isize, max_isize, max_samples=None):
        # initialize
        d = InsertSizeDistribution()
        d.min_isize = min_isize
        d.max_isize = max_isize
        d.arr = array.array('L', (0 for x in xrange(min_isize, max_isize+1)))
        count = 0
        for isize in sample_fragment_sizes(bamfh, genes, min_isize, max_isize):
            if (min_isize <= isize <= max_isize):
                # store in array
                d.arr[isize - min_isize] += 1
                count += 1
                if (max_samples is not None) and (count > max_samples):
                    break
        return d
