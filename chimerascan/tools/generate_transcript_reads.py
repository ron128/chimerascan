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

BASES_PER_LINE = 50

def split_seq(seq, chars_per_line):
    pos = 0
    newseq = []
    while pos < len(seq):
        if pos + chars_per_line > len(seq):        
            endpos = len(seq)
        else:
            endpos = pos + chars_per_line
        newseq.append(seq[pos:endpos])
        pos = endpos
    return '\n'.join(newseq)

def bed12_to_fasta(gene_feature_file, reference_seq_file):
    ref_fa = pysam.Fastafile(reference_seq_file)
    for g in GeneFeature.parse(open(gene_feature_file)):
        exon_seqs = []
        error_occurred = False
        for start, end in g.exons:
            seq = ref_fa.fetch(g.chrom, start, end)
            if not seq:
                error_occurred = True
                break
            exon_seqs.append(seq)
        if error_occurred:
            continue
        # make fasta record
        seq = ''.join(exon_seqs)
        if g.strand == '-':
            seq = DNA_reverse_complement(seq)
        # break seq onto multiple lines
        seqlines = split_seq(seq, BASES_PER_LINE)    
        yield (">%s range=%s:%d-%d gene=%s strand=%s\n%s" % 
               (config.GENE_REF_PREFIX + g.tx_name, g.chrom, start, end, g.gene_name, g.strand, seqlines))
    ref_fa.close()

def main():
    from optparse import OptionParser
    parser = OptionParser("usage: %prog [options] <index> <gene_features.txt> <outprefix>")
    parser.add_option("--error-rate", dest="error_rate", type="float", default=0.020)
    parser.add_option("-N", dest="num_reads", type="int", 
                      default=100, metavar="N", 
                      help="number of reads [default=%default]")
    parser.add_option("--rlen", dest="rlen", type="int", 
                      default=50, metavar="N", 
                      help="read length [default=%default]")
    parser.add_option("--isize", dest="isize", type="int", 
                      default=200, metavar="N",
                      help="insert size [default=%default]")
    parser.add_option("--isize-stdev", dest="isize_stdev", type="float",
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
    gene_feature_file = args[1]
    out_prefix = args[2]
    ref_fasta_file = os.path.join(index_dir, config.ALIGN_INDEX_FASTA_FILE)

    # make FASTA from the gene features in the input
    f = open("bubba", "w")
    for fasta_seq in bed12_to_fasta(gene_feature_file, ref_fasta_file):
        print >>f, fasta_seq
    f.close()
    # generate reads with 'wgsim'
    args = [os.path.join(options.wgsim_dir, "wgsim"),
            "-e", options.error_rate,
            "-d", options.isize,
            "-s", options.isize_stdev,
            "-N", options.num_reads,
            "-1", options.rlen,
            "-2", options.rlen,
            "-r", 0.0,
            "-R", 0.0,
            "bubba",
            out_prefix + "_1.fq",
            out_prefix + "_2.fq"]
    subprocess.call(map(str, args))
    os.remove("bubba")
                          
if __name__ == '__main__':
    main()
