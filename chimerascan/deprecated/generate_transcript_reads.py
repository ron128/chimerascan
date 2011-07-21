'''
Created on Jul 21, 2011

@author: mkiyer
'''
'''
Created on Jul 14, 2011

@author: mkiyer
'''
import os
import subprocess

import chimerascan.lib.config as config
from chimerascan.lib.feature import GeneFeature
from chimerascan.lib.seq import DNA_reverse_complement
from chimerascan import pysam

def main():
    from optparse import OptionParser
    parser = OptionParser("usage: %prog [options] <index> <outprefix>")
    parser.add_option("--error-rate", dest="error_rate", type="float", default=0.020)
    parser.add_option("-N", dest="num_reads", type="int", 
                      default=100, metavar="N", 
                      help="number of reads [default=%default]")
    parser.add_option("--rlen", dest="read_length", type="int", 
                      default=50, metavar="N", 
                      help="read length [default=%default]")
    parser.add_option("--isize", dest="insert_size", type="int", 
                      default=200, metavar="N",
                      help="insert size [default=%default]")
    parser.add_option("--isize-stdev", dest="insert_size_stdev", type="float",
                      default=20, metavar="N",
                      help="insert size standard deviation [defaul=%default]")
    parser.add_option("--library", dest="library_type", 
                      default="fr-unstranded",
                      help="library type [default=%default]")
    parser.add_option("--wgsim-dir", dest="wgsim_dir", 
                      default="", help="directory containing 'wgsim' tool "
                      "packaged with samtools [default=%default]")
    options, args = parser.parse_args()
    if len(args) < 2:
        parser.error("Not enough input arguments")
    # extract command line arguments
    index_dir = args[0]
    out_prefix = args[1]
    gene_feature_file = os.path.join(index_dir, config.GENE_FEATURE_FILE) 
    ref_fasta_file = os.path.join(index_dir, config.ALIGN_INDEX_FASTA_FILE)
    # generate reads with 'wgsim'
    args = [os.path.join(options.wgsim_path, "wgsim"),
            "-e", options.error_rate,
            "-d", options.isize,
            "-s", options.isize_stdev,
            "-N", options.num_reads,
            "-1", options.rlen,
            "-2", options.rlen,
            "-r", 0.0,
            "-R", 0.0,
            ref_fasta_file,
            out_prefix + "_1.fq",
            out_prefix + "_2.fq"]
            
            


    
    fh1 = open(out_prefix + "_1.fq", "w")
    fh2 = open(out_prefix + "_2.fq", "w")
    # iterate through genes
    fastafh = pysam.Fastafile(ref_fasta_file)
    for g in GeneFeature.parse(open(gene_feature_file)):
        rname = config.GENE_REF_PREFIX + g.tx_name
        seq = fastafh.fetch(reference=rname)
        # walk through seq and generate reads
        for end in xrange(insert_size, len(seq)):
            start = end - insert_size            
            seq1 = seq[start:(start+read_length)]
            seq2 = DNA_reverse_complement(seq[(end-read_length):end])
            tag = '%s:%d-%d' % (rname,start,end)
            print >>fh1, '@%s/1\n%s\n+\n%s' % (tag, seq1, 'I' * read_length)
            print >>fh2, '@%s/2\n%s\n+\n%s' % (tag, seq2, 'I' * read_length)    
    fastafh.close()
    fh1.close()
    fh2.close()
                                              
if __name__ == '__main__':
    main()
