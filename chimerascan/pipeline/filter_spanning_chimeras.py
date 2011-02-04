'''
Created on Jan 31, 2011

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
import logging
import operator
import subprocess
import tempfile
import os

# local lib imports
from chimerascan.lib import config
from chimerascan.lib.gene_to_genome2 import build_gene_to_genome_map, gene_to_genome_pos
from chimerascan.lib.stats import kl_divergence
# local imports
from nominate_chimeras import Chimera, MULTIMAP_BINS
        
class SpanningChimera(Chimera):
    def __init__(self):
        Chimera.__init__(self)
        self.spanning_reads = 0
        self.encomp_and_spanning = 0
        self.total_reads = 0
        self.junction_hist = None
        self.spanning_seq = None
        self.spanning_ids = None
        
    def from_list(self, fields):
        FIRST_COL = Chimera.LAST_COL + 1
        # get the chimera fields
        Chimera.from_list(self, fields)
        self.spanning_reads = int(fields[FIRST_COL])
        self.encomp_and_spanning = int(fields[FIRST_COL+1])
        self.total_reads = int(fields[FIRST_COL+2])
        self.junction_hist = map(int, fields[FIRST_COL+3].split(','))
        self.spanning_seq = fields[FIRST_COL+4]
        self.spanning_ids = fields[FIRST_COL+5]

    def to_list(self):
        fields = Chimera.to_list(self)
        fields.extend([self.spanning_reads, 
                       self.encomp_and_spanning, 
                       self.total_reads, 
                       ','.join(map(str, self.junction_hist)), 
                       self.spanning_seq, 
                       self.spanning_ids])
        return fields

    @staticmethod
    def parse(line_iter):
        for line in line_iter:
            fields = line.strip().split('\t')
            c = SpanningChimera()
            c.from_list(fields)
            yield c

def filter_multimapping(c, max_multimap=1, 
                        multimap_cov_ratio=0.0):
    '''
    generator that returns chimeras that based on the uniqueness of
    supporting reads.  chimeras supporting multimapping reads with more 
    than 'max_multimap' hits will be ignored, and chimeras with less
    than 'weighted_cov_ratio' fraction of coverage to reads will be ignored.

    for example, if a chimera has a coverage of 2.0, but has 200 reads,
    the ratio will be 2.0/200 = 1/100.  this suggests that the majority of
    reads supporting the candidate are multimapping.
    
    however, if there is one completely unique read supporting the candidate,
    then 1.0 out of 2.0 coverage is accountable to a single read.  so this
    candidate would pass the 'max_multimap' filter and not be removed
    '''
    # get index of first read
    for ind,x in enumerate(c.multimap_cov_hist):
        if x > 0:
            break
    mmap = MULTIMAP_BINS[ind]
    ratio = c.weighted_cov / float(c.encompassing_reads)
    if (mmap > max_multimap) or (ratio < multimap_cov_ratio):
        #logging.debug("Excluding chimera with %f cov, %d reads, and %s mmap hist" %
        #              (c.weighted_cov, c.encompassing_reads, c.multimap_cov_hist))
        return False
    return True

def filter_insert_size(c, max_isize):
    '''
    estimate the insert size by comparing the reads that map to the
    hypothetical 5'/3' transcript to the insert size distribution and
    remove chimeras that fail to meet this constraint
    '''
    if (c.mate5p.isize + c.mate3p.isize) <= (2*max_isize):
        return True
    else:
        #logging.warning("Removed %s due to insert size %d + %d > %d" %
        #                (c.name, c.mate5p.isize, c.mate3p.isize, 2*max_isize))
        return False

def filter_overlapping(c):
    return c.distance != 0

def build_junc_coverage_map(chimeras, ggmap):
    junc_cov_map = collections.defaultdict(lambda: [None, None, None])
    num_chimeras = 0
    for c in chimeras:
        num_chimeras += 1
        # convert to genomic coords
        # subtract one since 5' junc position is an open interval
        coord5p = gene_to_genome_pos(c.mate5p.tx_name, c.mate5p.end - 1, ggmap)
        coord3p = gene_to_genome_pos(c.mate3p.tx_name, c.mate3p.start, ggmap)
        # keep track of maximum coverage isoform
        pairkey = (coord5p, coord3p)
        paircov = (c.encomp_and_spanning, c.weighted_cov, c.total_reads)
        data = junc_cov_map[pairkey]
        if (data[0] is None) or (cmp(paircov, data[0]) > 0):
            # store encomp/spanning, weighted coverage, and total reads
            data[0] = paircov
            data[1] = c.mate5p.tx_name
            data[2] = c.mate3p.tx_name
    logging.debug("Parsed %d chimeras" % (num_chimeras))
    kept_isoforms = set(tuple(v[1:3]) for v in junc_cov_map.itervalues())
    #del junc_cov_map
    logging.debug("Kept %d highest coverage isoforms" % (len(kept_isoforms)))
    return kept_isoforms

def choose_highest_coverage_chimeras(input_file, ggmap):
    '''
    choose the highest coverage isoform pair using spanning reads,
    encompassing reads, and total reads as a measure.  ties will be
    broken by choosing a single gene pair arbitrarily 
    '''
    # break name into 5'/3' genes linked in a dictionary
    logging.debug("Building junction isoform coverage map")
    kept_isoforms_set = build_junc_coverage_map(SpanningChimera.parse(open(input_file)), ggmap)
    # write results
    logging.debug("Returning highest coverage chimeras")
    for c in SpanningChimera.parse(open(input_file)):
        pairkey = (c.mate5p.tx_name, c.mate3p.tx_name)
        if pairkey in kept_isoforms_set:
            yield c
    del kept_isoforms_set

def build_junc_permiscuity_map(chimeras, ggmap):
    junc5p_map = collections.defaultdict(lambda: collections.defaultdict(lambda: 0))
    junc3p_map = collections.defaultdict(lambda: collections.defaultdict(lambda: 0))
    for c in chimeras:
        # subtract one since 5' junc position is an open interval
        coord5p = gene_to_genome_pos(c.mate5p.tx_name, c.mate5p.end - 1, ggmap)
        coord3p = gene_to_genome_pos(c.mate3p.tx_name, c.mate3p.start, ggmap)
        # keep track of total reads eminating from each 5' junction
        # by keeping a dictionary for each 5' junction to all 3' junctions
        # that stores the maximum coverage at that 5'/3' pair
        partners = junc5p_map[coord5p]
        count = partners[coord3p]
        partners[coord3p] = max(count, c.weighted_cov)
        # repeat for 3' partner
        partners = junc3p_map[coord3p]
        count = partners[coord5p]
        partners[coord5p] = max(count, c.weighted_cov)
        #print '5P', c.mate5p.gene_name, len(partners), sum(partners.itervalues())
        #print '3P', c.mate3p.gene_name, len(partners), sum(partners.itervalues())
    return junc5p_map, junc3p_map

def collect_permiscuity_stats(input_file, ggmap):
    # break name into 5'/3' genes linked in a dictionary
    logging.debug("Building chimera permiscuity map")
    juncmap5p, juncmap3p = \
        build_junc_permiscuity_map(SpanningChimera.parse(open(input_file)), ggmap)
    return juncmap5p, juncmap3p

def calc_permiscuity(c, juncmap5p, juncmap3p, ggmap):
    # subtract one since 5' junc position is an open interval
    coord5p = gene_to_genome_pos(c.mate5p.tx_name, c.mate5p.end - 1, ggmap)
    coord3p = gene_to_genome_pos(c.mate3p.tx_name, c.mate3p.start, ggmap)
    partners = juncmap5p[coord5p]
    cov = partners[coord3p]
    total_cov = sum(partners.itervalues())
    frac5p = cov / float(total_cov)
    partners = juncmap3p[coord3p]
    cov = partners[coord5p]
    total_cov = sum(partners.itervalues())
    frac3p = cov / float(total_cov)
    return frac5p, frac3p

def make_temp(base_dir, suffix=''):
    fd,name = tempfile.mkstemp(suffix=suffix, prefix='tmp', dir=base_dir)
    os.close(fd)
    return name

def filter_spanning_chimeras(input_file, output_file, gene_file):
#    tmpfile1 = make_temp(base_dir=os.path.dirname(output_file),
#                         suffix='.bedpe')
#    logging.debug("Filtering chimeras")
#    fh = open(tmpfile1, "w")
#    for c in SpanningChimera.parse(open(input_file)):
#        res = filter_multimapping(c, max_multimap=max_multimap,
#                                  multimap_cov_ratio=multimap_cov_ratio)
#        res = res and filter_insert_size(c, max_isize)
#        res = res and filter_overlapping(c)
#        if not res:
#            continue
#        print >>fh, '\t'.join(map(str, c.to_list()))
#    fh.close()
    logging.debug("Building gene/genome index")
    ggmap = build_gene_to_genome_map(open(gene_file))
    logging.debug("Choosing highest coverage chimeras")
    tmpfile1 = make_temp(base_dir=os.path.dirname(output_file),
                         suffix='.bedpe')
    fh = open(tmpfile1, "w")
    for c in choose_highest_coverage_chimeras(input_file, ggmap):
        print >>fh, '\t'.join(map(str, c.to_list()))
    fh.close()
    logging.debug("Finding junction permiscuity")
    juncmap5p, juncmap3p = collect_permiscuity_stats(tmpfile1, ggmap)
    fh = open(output_file, "w")
    for c in SpanningChimera.parse(open(tmpfile1)):
        frac5p, frac3p = calc_permiscuity(c, juncmap5p, juncmap3p, ggmap)
        kldiv = kl_divergence(c.junction_hist)
        print >>fh, '\t'.join(['\t'.join(map(str,c.to_list())), str(kldiv), 
                               str(frac5p), str(frac3p)])
    fh.close()
    # delete tmp files
    os.remove(tmpfile1)


def main():
    from optparse import OptionParser
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = OptionParser("usage: %prog [options] <sortedchimeras.bedpe> <chimeras.txt>")
    parser.add_option("--index", dest="index_dir",
                      help="Path to chimerascan index directory")
    options, args = parser.parse_args()
    gene_file = os.path.join(options.index_dir, config.GENE_FEATURE_FILE)
    input_file = args[0]
    output_file = args[1]
    filter_spanning_chimeras(input_file, output_file, gene_file)

if __name__ == "__main__":
    main()