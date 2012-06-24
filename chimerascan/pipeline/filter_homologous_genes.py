'''
Created on Aug 1, 2011

@author: mkiyer
'''
import logging
import os
import collections
import subprocess

import pysam

from chimerascan.lib import config
from chimerascan.lib.chimera import Chimera
from chimerascan.bx.intersection import IntervalTree, Interval

def get_mapped_read_intervals(c, min_isize, max_isize, homolog_segment_length):
    start5p = max(0, c.tx_end_5p - min_isize + homolog_segment_length)         
    end5p = max(0, c.tx_end_5p + max_isize - homolog_segment_length)
    if start5p > end5p:
        end5p = start5p + homolog_segment_length    
    start3p = max(0, c.tx_start_3p - max_isize + homolog_segment_length)
    end3p = max(0, c.tx_start_3p + min_isize - homolog_segment_length)  
    if start3p > end3p:
        end3p = start3p + homolog_segment_length
    return start5p, end5p, start3p, end3p
    
def filter_homologous_genes(input_file, 
                            output_file, 
                            index_dir,
                            homolog_segment_length,
                            min_isize,
                            max_isize,
                            maxhits,
                            num_processors,
                            tmp_dir):
    logging.debug("Parameters")
    logging.debug("\thomolog segment length: %d" % (homolog_segment_length))
    logging.debug("\tmin fragment size: %d" % (min_isize))
    logging.debug("\tmax fragment size: %d" % (max_isize))
    # open the reference sequence fasta file
    ref_fasta_file = os.path.join(index_dir, config.TRANSCRIPTOME_INDEX + ".fa")
    ref_fa = pysam.Fastafile(ref_fasta_file)
    interval_trees_3p = collections.defaultdict(lambda: IntervalTree())
    # generate FASTA file of sequences to use in mapping
    logging.debug("Generating homologous sequences to test")
    fasta5p = os.path.join(tmp_dir, "homologous_5p.fa")    
    f = open(fasta5p, "w")
    for c in Chimera.parse(open(input_file)):
        start5p, end5p, start3p, end3p = get_mapped_read_intervals(c, min_isize, max_isize, homolog_segment_length)
        # add 3' gene to interval trees
        interval_trees_3p[c.tx_name_3p].insert_interval(Interval(start3p, end3p, value=c.name))
        # extract sequence of 5' gene
        seq5p = ref_fa.fetch(c.tx_name_5p, start5p, end5p)
        for i in xrange(0, len(seq5p) - homolog_segment_length):
            print >>f, ">%s,%s:%d-%d\n%s" % (c.name,c.tx_name_5p,
                                             start5p+i,
                                             start5p+i+homolog_segment_length,
                                             seq5p[i:i+homolog_segment_length])
    f.close()
    # map 5' sequences to reference using bowtie
    logging.debug("Mapping homologous sequences")
    bowtie2_index = os.path.join(index_dir, config.TRANSCRIPTOME_INDEX)
    sam5p = os.path.join(tmp_dir, "homologous_5p.sam")
    args = [config.BOWTIE2_BIN, 
            '-p', num_processors, '--phred33',
            '--end-to-end', '--very-sensitive', '--reorder',
            '-f', '-k', maxhits,
            '-x', bowtie2_index,
            '-U', fasta5p,
            "-S", sam5p]
    retcode = subprocess.call(map(str,args))
    if retcode != 0:
        return config.JOB_ERROR
    # analyze results for homologous genes
    logging.debug("Analyzing mapping results")
    samfh = pysam.Samfile(sam5p, "r")
    tid_rname_map = dict((i,refname) for i,refname in enumerate(samfh.references))
    homologous_chimeras = set()
    for r in pysam.Samfile(sam5p, "r"):
        if r.is_unmapped:
            continue
        # reference name must be in list of 3' chimeras
        rname = tid_rname_map[r.tid]        
        if rname not in interval_trees_3p:
            continue
        # get chimera name from 'qname'
        chimera_name = r.qname.split(",")[0]
        for hit in interval_trees_3p[rname].find(r.pos,r.aend):
            if hit.value == chimera_name:
                homologous_chimeras.add(chimera_name)
    # write output
    logging.debug("Writing output")
    f = open(output_file, "w")
    for c in Chimera.parse(open(input_file)):
        if c.name in homologous_chimeras:
            logging.debug("Removing homologous chimera %s between %s and %s" % 
                          (c.name, c.gene_name_5p, c.gene_name_3p))
            continue
        print >>f, '\t'.join(map(str, c.to_list()))        
    f.close()
    # cleanup
    if os.path.exists(fasta5p):
        os.remove(fasta5p)
    if os.path.exists(sam5p):
        os.remove(sam5p)    
    return config.JOB_SUCCESS
    

def main():
    from optparse import OptionParser
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = OptionParser("usage: %prog [options] <index_dir> "
                          "<in.txt> <out.txt>")
    parser.add_option("--homolog-segment-length", dest="homolog_segment_length",
                      type="int", default=25, 
                      help="Segment length to consider when searching for "
                      "homologous regions [default=%default]")
    parser.add_option('--min-fragment-length', dest="min_fragment_length", 
                      type="int", default=100)
    parser.add_option('--max-fragment-length', dest="max_fragment_length", 
                      type="int", default=300)    
    parser.add_option("--bowtie-bin", dest="bowtie_bin",
                      default="bowtie", 
                      help="Path to bowtie binary [default: %default]")
    parser.add_option("-p", type="int", dest="num_processors", default=1,
                      help="Number of processors to use [default: %default]")
    parser.add_option("--tmp-dir", dest="tmp_dir",
                      default=".", 
                      help="Temporary directory [default=%default]")
    options, args = parser.parse_args()
    index_dir = args[0]
    input_file = args[1]
    output_file = args[2]
    return filter_homologous_genes(input_file, output_file, index_dir,
                                   homolog_segment_length=options.homolog_segment_length,
                                   min_isize=options.min_fragment_length,
                                   max_isize=options.max_fragment_length,
                                   bowtie_bin=options.bowtie_bin,
                                   num_processors=options.num_processors,
                                   tmp_dir=options.tmp_dir)

if __name__ == "__main__":
    main()