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
import os

# local lib imports
from chimerascan.lib import config
from chimerascan.lib.gene_to_genome2 import build_gene_to_genome_map, gene_to_genome_pos
from chimerascan.lib.stats import binomial_cdf

# local imports
from merge_spanning_alignments import SpanningChimera

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
        paircov = (c.encomp_and_spanning, c.weighted_cov, c.encomp_or_spanning)
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

def filter_read_balance(c, pval):
    '''
    returns True if binomial test pvalue for strand balance is greater than
    'pval', False otherwise
    '''    
    # count reads on either strand
    mate_counts = [0, 0]
    for r in c.spanning_reads:
        mate_counts[r.mate] += 1
    p = binomial_cdf(0.5, c.num_spanning_reads, min(mate_counts))        
    if p <= pval:
        logging.warning("Filtered chimera spanning reads=%d read1=%d read2=%d pval=%f" %
                        (c.num_spanning_reads, mate_counts[0], mate_counts[1], p))
    return p > pval

def filter_insert_size(c, max_isize):
    '''
    estimate the insert size by comparing the reads that map to the
    hypothetical 5'/3' transcript to the insert size distribution and
    remove chimeras that fail to meet this constraint

    returns True if chimera agrees with insert size distribution, 
    false otherwise
    '''
    if max_isize <= 0: 
        return True
    if (c.mate5p.isize > max_isize) or (c.mate3p.isize > max_isize):
        logging.warning("Removed %s due to insert size %d + %d > %d" %
                        (c.name, c.mate5p.isize, c.mate3p.isize, 2*max_isize))
        return False
    return True

def filter_spanning_chimeras(input_file, output_file, gene_file,
                             mate_pval, max_isize):
    '''
    processes chimera isoforms and chooses the one with the 
    highest coverage and omits the rest
    '''
    logging.debug("Building gene/genome index")
    ggmap = build_gene_to_genome_map(open(gene_file))
    logging.debug("Choosing highest coverage chimeras")
    fh = open(output_file, "w")
    for c in choose_highest_coverage_chimeras(input_file, ggmap):
        res = True
        res = res and filter_insert_size(c, max_isize)
        if res:
            print >>fh, '\t'.join(['\t'.join(map(str,c.to_list()))])
    fh.close()

def main():
    from optparse import OptionParser
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = OptionParser("usage: %prog [options] <sortedchimeras.bedpe> <chimeras.txt>")
    parser.add_option("--index", dest="index_dir",
                      help="Path to chimerascan index directory")
    parser.add_option("--mate-pval", type="float", metavar="p", 
                      dest="mate_pval", default=0.01,                       
                      help="p-value to reject chimera based on binomial "
                      "test that balance of read1/read2 "
                      "should be 50/50 [default=%default]")
    parser.add_option("--max-isize", type="float", dest="max_isize",
                      default=-1, help="Maximum predicted insert size of "
                      "fragments spanning a hypothetical chimeric junction "
                      "[default=%default]")    
    options, args = parser.parse_args()
    gene_file = os.path.join(options.index_dir, config.GENE_FEATURE_FILE)
    input_file = args[0]
    output_file = args[1]
    filter_spanning_chimeras(input_file, output_file, gene_file,
                             mate_pval=options.mate_pval,
                             max_isize=options.max_isize)

if __name__ == "__main__":
    main()