'''
Created on Feb 24, 2012

@author: mkiyer
'''
import logging
import os
import sys

from chimerascan import pysam
from chimerascan.lib import config
from chimerascan.lib.base import LibraryTypes

def extract_tophat_encompassing_reads(index_dir, tophat_bam_file, 
                                      encompassing_bam_file,
                                      max_isize, library_type):
    gene_file = os.path.join(index_dir, config.GENE_FEATURE_FILE)
    bamfh = pysam.Samfile(tophat_bam_file, "rb")
    for r in bamfh:
        if (r.is_unmapped) or (r.mate_is_unmapped):
            continue
        if r.rname != r.mrnm:
            print r.qname, r.rname, r.pos, r.is_reverse, r.mrnm, r.mpos, r.mate_is_reverse       
    bamfh.close()
    

def main():
    from optparse import OptionParser
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = OptionParser("usage: %prog [options] <index> <tophat.bam> <tophat_encompassing.bam>")
    parser.add_option('--max-fragment-length', dest="max_fragment_length", 
                      type="int", default=1000)
    parser.add_option('--library', dest="library_type", 
                      default=LibraryTypes.FR_UNSTRANDED)
    options, args = parser.parse_args()    
    index_dir = args[0]
    tophat_bam_file = args[1]
    encompassing_bam_file = args[2]
    return extract_tophat_encompassing_reads(index_dir, tophat_bam_file, 
                                             encompassing_bam_file,
                                             max_isize=options.max_fragment_length,
                                             library_type=options.library_type)

if __name__ == '__main__':
    sys.exit(main())