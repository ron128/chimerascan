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
import logging
import os

def filter_spurious_spanning_reads(c):
    """
    removes spanning reads that lack paired 
    """
    pass

def filter_coverage(c, encomp_min=2, spanning_threshold=1):
    """
    filters chimeras with fewer than 'spanning_threshold' spanning
    reads and fewer than 'encomp_min' encompassing reads
    """
    
    
    pass


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
                      "[default=%default]")
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