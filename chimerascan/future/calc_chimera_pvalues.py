'''
Created on Aug 2, 2011

@author: mkiyer
'''
import logging
import os
import re

from chimerascan import pysam
from chimerascan.lib import config
from chimerascan.lib.chimera import Chimera
from chimerascan.lib.batch_sort import batch_sort
from chimerascan.lib.fragment_size_distribution import InsertSizeDistribution

from chimerascan.pipeline.resolve_discordant_reads import make_discordant_read_stats_file, sort_read_stats_by_read_name, group_by_attr, ChimeraStats

_bowtie_num_reads_re = r'^# reads with at least one reported alignment:\s+(\d+)'

def get_bowtie_num_mapped_reads(log_file):
    for line in open(log_file):
        line = line.strip()
        m = re.match(_bowtie_num_reads_re, line)
        if m is not None:
            return int(m.group(1))
    return None

def calc_percent_discordant_reads(input_file, bowtie_log_file, isize_dist, min_isize_prob, tmp_dir):
    #
    # count total number of mapped reads
    #
    num_mapped_frags = get_bowtie_num_mapped_reads(bowtie_log_file)
    if num_mapped_frags is None:
        return config.JOB_ERROR
    logging.debug("\tmapped reads: %d" % (num_mapped_frags))
    #
    # parse chimeras and output read statistics to a file
    #
    logging.debug("Getting discordant read information")
    read_stats_file = os.path.join(tmp_dir, "read_stats.txt")
    make_discordant_read_stats_file(input_file, read_stats_file, isize_dist)
    #
    # now sort the read/chimera stats list
    #
    logging.debug("Sorting reads by read name")
    sorted_read_stats_file = os.path.join(tmp_dir, "read_stats.rname_sorted.txt")
    sort_read_stats_by_read_name(read_stats_file, sorted_read_stats_file, tmp_dir)
    #
    # count number of unique discordant reads that map within insert size 
    # probability threshold
    #
    logging.debug("Counting valid discordant reads")
    num_discordant_reads = 0
    num_discordant_reads_within_isize_range = 0
    for rname,readstats in group_by_attr(ChimeraStats.parse(open(sorted_read_stats_file)), 'qname'):
        num_discordant_reads += 1        
        max_isize_prob = max(s.isize_prob for s in readstats)
        if max_isize_prob >= min_isize_prob:
            num_discordant_reads_within_isize_range += 1
    logging.debug("\tdiscordant reads: %d" % (num_discordant_reads))
    logging.debug("\tdiscordant reads within insert size range: %d" % 
                  (num_discordant_reads_within_isize_range))
    # remove temporary files
    os.remove(read_stats_file)
    os.remove(sorted_read_stats_file)    
    return (num_mapped_frags + num_discordant_reads), num_discordant_reads_within_isize_range

def calc_chimera_pvalues(input_file,
                         bam_file, 
                         num_mapped_reads, 
                         num_discordant_reads_within_isize_range):
    # calc discordant reads per million
    percent_discordant = num_discordant_reads_within_isize_range / float(num_mapped_reads)
    # open BAM file for checking wild-type isoforms
    bamfh = pysam.Samfile(bam_file, "rb")
    for c in Chimera.parse(open(input_file)):        
        # count 5' and 3' reads
        rname5p = config.GENE_REF_PREFIX + c.tx_name_5p
        rname3p = config.GENE_REF_PREFIX + c.tx_name_3p        
        num_reads_5p = len(set(r.qname for r in bamfh.fetch(rname5p, c.tx_start_5p, c.tx_end_5p)))
        num_reads_3p = len(set(r.qname for r in bamfh.fetch(rname3p, c.tx_start_3p, c.tx_end_3p)))
        # expected number of discordant reads
        exp_discordant_5p = num_reads_5p * percent_discordant
        exp_discordant_3p = num_reads_3p * percent_discordant
        print c.gene_name_5p, c.gene_name_3p, num_reads_5p, num_reads_3p, exp_discordant_5p, exp_discordant_3p
    bamfh.close()    

def main():
    import sys
    calc_chimera_pvalues(sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]))
    return
    from optparse import OptionParser
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = OptionParser("usage: %prog [options] <in.txt> <bowtie.log> <isizedist.txt>")
    parser.add_option("--min-isize-prob", dest="min_isize_prob", 
                      type="float", default=0.01)
    options, args = parser.parse_args()
    input_file = args[0]
    bowtie_log_file = args[1]
    isize_dist_file = args[2]
    # read insert size distribution
    isize_dist = InsertSizeDistribution.from_file(open(isize_dist_file))
    calc_percent_discordant_reads(input_file, bowtie_log_file, isize_dist,
                                  min_isize_prob=options.min_isize_prob, 
                                  tmp_dir=".")

if __name__ == '__main__':
    main()