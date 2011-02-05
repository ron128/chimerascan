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
import tempfile
import os

# local lib imports
from chimerascan.lib import config
from chimerascan.lib.gene_to_genome2 import build_gene_to_genome_map, gene_to_genome_pos
from chimerascan.lib.stats import binomial_cdf
# local imports
from nominate_chimeras import Chimera, MULTIMAP_BINS
        
def filter_multimapping(c, max_multimap=1, multimap_cov_ratio=0.0):
    '''
    returns True/False based on the uniqueness of supporting reads.  
    chimeras with multimapping reads with more than 'max_multimap' 
    hits will be ignored, and chimeras with less than 'weighted_cov_ratio' 
    fraction of coverage to reads will be ignored.

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
    if (mmap <= max_multimap) and (ratio >= multimap_cov_ratio):
        return True
    if c.weighted_cov >= 5:
        logging.debug("Excluding chimera with %f cov, %d reads, and %s mmap hist" %
                      (c.weighted_cov, c.encompassing_reads, c.multimap_cov_hist))
    return False

def filter_insert_size(c, max_isize):
    '''
    estimate the insert size by comparing the reads that map to the
    hypothetical 5'/3' transcript to the insert size distribution and
    remove chimeras that fail to meet this constraint

    returns True if chimera agrees with insert size distribution, 
    false otherwise
    '''
    if (c.mate5p.isize + c.mate3p.isize) <= (2*max_isize):
        return True
    #logging.warning("Removed %s due to insert size %d + %d > %d" %
    #                (c.name, c.mate5p.isize, c.mate3p.isize, 2*max_isize))
    return False

def filter_overlapping(c):
    '''
    filter chimeras on overlapping genes
    
    returns True if chimera is not overlapping, False otherwise
    '''
    return c.distance != 0

def filter_strand_balance(c, pval):
    '''
    returns True if binomial test pvalue for strand balance is greater than
    'pval', False otherwise
    '''
    p = binomial_cdf(0.5, c.encompassing_reads, min(c.strand_reads))        
    if p <= pval:
        logging.warning("Filtered chimera reads=%d '+'=%d '-'=%d pval=%f" %
                        (c.encompassing_reads, c.strand_reads[0], 
                         c.strand_reads[1], p))
    return p > pval

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
        build_junc_permiscuity_map(Chimera.parse(open(input_file)), ggmap)
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

def filter_encompassing_chimeras(input_file, output_file, gene_file,
                                 max_multimap=1,
                                 multimap_cov_ratio=0.10,
                                 max_isize=1000,
                                 strand_pval=0.01,
                                 keep_overlap=False):
    logging.debug("Filtering chimeras")
    logging.debug("Must have a read with <= %d multimaps" % (max_multimap))
    logging.debug("Coverage to reads ratio >= %f" % (multimap_cov_ratio))
    logging.debug("Insert size < %d" % (max_isize))
    logging.debug("Strand balance p-value > %f" % (strand_pval))
    # first perform basic filtering
    tmpfile1 = make_temp(base_dir=os.path.dirname(output_file),
                         suffix='.bedpe')
    fh = open(tmpfile1, "w")
    for c in Chimera.parse(open(input_file)):
        res = filter_multimapping(c, max_multimap=max_multimap, 
                                  multimap_cov_ratio=multimap_cov_ratio)
        res = res and filter_insert_size(c, max_isize)
        if not keep_overlap:
            res = res and filter_overlapping(c)
        res = res and filter_strand_balance(c, strand_pval)
        if res:
            print >>fh, '\t'.join(map(str, c.to_list()))
    fh.close()
    logging.debug("Building gene/genome index")
    ggmap = build_gene_to_genome_map(open(gene_file))
    logging.debug("Finding junction permiscuity")
    juncmap5p, juncmap3p = collect_permiscuity_stats(tmpfile1, ggmap)
    fh = open(output_file, "w")
    for c in Chimera.parse(open(tmpfile1)):
        frac5p, frac3p = calc_permiscuity(c, juncmap5p, juncmap3p, ggmap)
        c.mate5p.frac = frac5p
        c.mate3p.frac = frac3p
        print >>fh, '\t'.join(map(str, c.to_list()))
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
    parser.add_option("--max-multimap", type="int", dest="max_multimap", 
                      default=1, help="Threshold to eliminate multimapping "
                      "chimeras, where '1' is completely unique, '2' is "
                      "multimapping to two locations, etc.")
    parser.add_option("--multimap-ratio", type="float", dest="multimap_cov_ratio",
                      default=0.10, help="Ratio of weighted coverage to "
                      "total encompassing reads below which chimeras are "
                      "considered false positives and removed "
                      "[default=%default]")
    parser.add_option("--max-isize", type="float", dest="max_isize",
                      default=500, help="Maximum predicted insert size of "
                      "fragments spanning a hypothetical chimeric junction "
                      "[default=%default]")
    parser.add_option("--strand-pval", type="float", metavar="p", 
                      dest="strand_pval", default=0.01,                       
                      help="p-value to reject chimera based on binomial "
                      "test that balance of +/- strand encompassing reads "
                      "should be 50/50 [default=%default]")
    parser.add_option("--keep-overlap", action="store_true", 
                      default=False, dest="keep_overlap",
                      help="keep chimera candidates that occur between "
                      "overlapping genes.  these are likely to be splice "
                      "variants that did not occur in the reference. "
                      "[default=%default")
    options, args = parser.parse_args()
    gene_file = os.path.join(options.index_dir, config.GENE_FEATURE_FILE)
    input_file = args[0]
    output_file = args[1]
    filter_encompassing_chimeras(input_file, output_file, gene_file,
                                 max_multimap=options.max_multimap,
                                 multimap_cov_ratio=options.multimap_cov_ratio,
                                 max_isize=options.max_isize,
                                 strand_pval=options.strand_pval,
                                 keep_overlap=options.keep_overlap)

if __name__ == "__main__":
    main()