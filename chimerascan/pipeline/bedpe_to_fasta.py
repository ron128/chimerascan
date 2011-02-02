'''
Created on Oct 26, 2010

@author: mkiyer
'''
import logging
import collections
import os

# local imports
from chimerascan import pysam
from chimerascan.lib import config

def bedpe_to_junction_fasta(bed_file, reference_seq_file, read_length,
                            fasta_output_fh, junc_output_fh):
    gene_fasta_prefix = config.GENE_REF_PREFIX
    ref_fa = pysam.Fastafile(reference_seq_file)
    juncs = collections.defaultdict(lambda: [])
    for line in open(bed_file):
        #print line
        fields = line.strip().split('\t')
        ref5p, start5p, end5p = fields[0], int(fields[1]), int(fields[2])
        ref3p, start3p, end3p = fields[3], int(fields[4]), int(fields[5])
        # join end of 5' ref with beginning of 3' ref
        junc_start5p = max(start5p, end5p - read_length + 1)
        junc_end3p = min(end3p, start3p + read_length - 1)
        # fetch sequence
        seq5p = ref_fa.fetch(gene_fasta_prefix + ref5p, junc_start5p, end5p)
        seq3p = ref_fa.fetch(gene_fasta_prefix + ref3p, start3p, junc_end3p)
        seq = seq5p + seq3p
        if len(seq) < (read_length*2) - 2:
            logging.warning("Could not extract sequence of length >%d from BEDPE, only retrieved sequence of (%d,%d) for gene %s" % 
                            ((read_length*2)-2, len(seq5p), len(seq3p), line.strip()))
        # add sequence to dictionary and group fusion candidates together
        # if they have the same junction sequence
        juncs[(len(seq5p),seq)].append(fields)
    # now extract the unique junction sequences
    # and write them to a fasta file
    junc_index = 1    
    for junc_seq_info,bedpe_fields in juncs.iteritems():
        left_seq_length, seq = junc_seq_info
        junc_name = "BLABBY%07d:%d" % (junc_index,left_seq_length)
        print >>fasta_output_fh, ">%s\n%s" % (junc_name, seq)
        for fields in bedpe_fields:
            print >>junc_output_fh, '\t'.join([junc_name] + fields)
        junc_index += 1

def main():
    from optparse import OptionParser
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = OptionParser("usage: %prog [options] <chimeras.bedpe> <out.fasta> <out.juncs>")
    parser.add_option("--rlen", type="int", dest="read_length", default=None)
    parser.add_option("--index", dest="index_dir") 
    options, args = parser.parse_args()
    bedpe_file = args[0]
    output_fasta_file = args[1]
    output_junc_file = args[2]
    ref_fasta_file = os.path.join(options.index_dir, config.ALIGN_INDEX + ".fa")    
    bedpe_to_junction_fasta(bedpe_file, ref_fasta_file,
                            options.read_length,
                            open(output_fasta_file, "w"),
                            open(output_junc_file, "w"))

if __name__ == '__main__':
    main()

