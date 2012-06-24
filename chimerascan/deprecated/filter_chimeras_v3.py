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
import collections
import os

import pysam

from chimerascan.lib.transcriptome import build_transcript_genome_map, \
    transcript_to_genome_pos, build_transcript_cluster_map
from chimerascan.lib.chimera import Chimera
from chimerascan.lib.feature import TranscriptFeature
from chimerascan.lib import config

def filter_unique_frags(c, threshold):
    """
    filters chimeras with less than 'threshold' unique
    alignment positions supporting the chimera 
    """
    return c.get_num_unique_positions() >= threshold

def get_wildtype_frags_5p(rname, start, end, bamfh):
    num_wildtype_frags = len(set(r.qname for r in bamfh.fetch(rname, start, end)
                                 if (not r.mate_is_unmapped) and (r.pnext >= end)))
    return num_wildtype_frags

def get_wildtype_frags_3p(rname, start, end, bamfh):
    num_wildtype_frags = len(set(r.qname for r in bamfh.fetch(rname, start, end)
                                 if (not r.mate_is_unmapped) and (r.pnext < start)))
    return num_wildtype_frags

def get_wildtype_frags(c, bamfh):
    num_wt_frags_5p = get_wildtype_frags_5p(c.tx_name_5p, c.tx_start_5p, c.tx_end_5p, bamfh)
    num_wt_frags_3p = get_wildtype_frags_3p(c.tx_name_3p, c.tx_start_3p, c.tx_end_3p, bamfh)
    return num_wt_frags_5p, num_wt_frags_3p

def filter_chimeric_isoform_fraction(c, frac, bamfh):
    """
    filters chimeras with fewer than 'threshold' total
    unique read alignments
    """
    num_wt_frags_5p, num_wt_frags_3p = get_wildtype_frags(c, bamfh)
    num_chimeric_frags = c.get_num_frags()
    ratio5p = float(num_chimeric_frags) / (num_chimeric_frags + num_wt_frags_5p)
    ratio3p = float(num_chimeric_frags) / (num_chimeric_frags + num_wt_frags_3p)
    #print c.gene_name_5p, c.gene_name_3p, "chimeras", num_chimeric_frags, "wt5p", num_wt_frags_5p, "wt3p", num_wt_frags_3p, "r5p", ratio5p, "r3p", ratio3p
    return min(ratio5p, ratio3p) >= frac

def read_false_pos_file(filename):
    false_pos_chimeras = set()
    for line in open(filename):
        fields = line.strip().split("\t")
        tx_name_5p, end5p, tx_name_3p, start3p = fields
        end5p = int(end5p)
        start3p = int(start3p)
        false_pos_chimeras.add((tx_name_5p, end5p, tx_name_3p, start3p))
    return false_pos_chimeras

def filter_encompassing_chimeras(input_file, output_file, min_frags):
    num_chimeras = 0
    num_filtered_chimeras = 0
    f = open(output_file, "w") 
    for c in Chimera.parse(open(input_file)):
        num_chimeras += 1
        if c.get_num_frags() < min_frags:
            continue
        num_filtered_chimeras += 1
        print >>f, '\t'.join(map(str, c.to_list()))
    f.close()
    logging.debug("\tchimeras: %d" % (num_chimeras))
    logging.debug("\tfiltered chimeras: %d" % (num_filtered_chimeras))
    return config.JOB_SUCCESS

def filter_chimeras(input_file, output_file,
                    index_dir, bam_file,
                    unique_frags,
                    isoform_fraction,
                    false_pos_file):
    logging.debug("Parameters")
    logging.debug("\tunique fragments: %f" % (unique_frags))
    logging.debug("\tfraction of wild-type isoform: %f" % (isoform_fraction))
    logging.debug("\tfalse positive chimeras file: %s" % (false_pos_file))
    # get false positive chimera list
    if (false_pos_file is not None) and (false_pos_file is not ""):
        logging.debug("Loading false positive chimeras")
        false_pos_pairs = read_false_pos_file(false_pos_file)
    else:
        false_pos_pairs = set()
    # open BAM file for checking wild-type isoform
    bamfh = pysam.Samfile(bam_file, "rb")
    # filter chimeras
    logging.debug("Filtering chimeras")
    num_chimeras = 0
    num_filtered_chimeras = 0    
    f = open(output_file, "w")   
    for c in Chimera.parse(open(input_file)):
        num_chimeras += 1
        good = filter_unique_frags(c, unique_frags)
        if not good:
            continue          
        false_pos_key = (c.tx_name_5p, c.tx_end_5p, 
                         c.tx_name_3p, c.tx_start_3p)
        good = good and (false_pos_key not in false_pos_pairs)
        if not good:
            continue
        good = good and filter_chimeric_isoform_fraction(c, isoform_fraction, bamfh)        
        if good:
            print >>f, '\t'.join(map(str, c.to_list()))
            num_filtered_chimeras += 1
    f.close()
    logging.debug("Total chimeras: %d" % num_chimeras)
    logging.debug("Filtered chimeras: %d" % num_filtered_chimeras)
    # cleanup memory for false positive chimeras
    del false_pos_pairs
    bamfh.close()
    return config.JOB_SUCCESS

def get_highest_coverage_isoforms(input_file, transcripts):
    # build lookup from transcript name to cluster id
    transcript_cluster_map = dict((str(t.tx_id),t.cluster_id) for t in transcripts)
    # build a lookup table to get genome coordinates from transcript 
    # coordinates
    transcript_genome_map = build_transcript_genome_map(transcripts)
    cluster_chimera_dict = collections.defaultdict(lambda: [])
    for c in Chimera.parse(open(input_file)):
        # TODO: adjust this to score chimeras differently!
        key = (c.name, c.get_num_frags())
        # get cluster of overlapping genes
        cluster5p = transcript_cluster_map[c.tx_name_5p]
        cluster3p = transcript_cluster_map[c.tx_name_3p]
        # get genomic positions of breakpoints
        coord5p = transcript_to_genome_pos(c.tx_name_5p, c.tx_end_5p-1, transcript_genome_map)
        coord3p = transcript_to_genome_pos(c.tx_name_3p, c.tx_start_3p, transcript_genome_map)
        # add to dictionary
        cluster_chimera_dict[(cluster5p,cluster3p,coord5p,coord3p)].append(key)    
    # choose highest coverage chimeras within each pair of clusters
    logging.debug("Finding highest coverage isoforms")
    kept_chimeras = set()
    for stats_list in cluster_chimera_dict.itervalues():
        stats_dict = collections.defaultdict(lambda: set())
        for stats_info in stats_list:
            # index chimera names
            stats_dict[stats_info[1:]].add(stats_info[0])
        # find highest scoring key
        sorted_keys = sorted(stats_dict.keys(), reverse=True)
        kept_chimeras.update(stats_dict[sorted_keys[0]])
    return kept_chimeras

def filter_highest_coverage_isoforms(index_dir, input_file, output_file):
    # read transcripts
    logging.debug("Reading transcripts")
    transcript_file = os.path.join(index_dir, config.TRANSCRIPT_FEATURE_FILE)
    transcripts = list(TranscriptFeature.parse(open(transcript_file)))
    # find highest coverage chimeras among isoforms
    kept_chimeras = get_highest_coverage_isoforms(input_file, transcripts)
    num_filtered_chimeras = 0
    f = open(output_file, "w")
    for c in Chimera.parse(open(input_file)):
        if c.name in kept_chimeras:
            num_filtered_chimeras += 1
            print >>f, '\t'.join(map(str, c.to_list()))
    f.close()
    logging.debug("\tAfter choosing best isoform: %d" % 
                  num_filtered_chimeras)
    return config.JOB_SUCCESS

def main():
    from optparse import OptionParser
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = OptionParser("usage: %prog [options] <index_dir> "
                          "<sorted_aligned_reads.bam> <in.txt> <out.txt>")
    parser.add_option("--unique-frags", type="float", default=2.0,
                      dest="unique_frags", metavar="N",
                      help="Filter chimeras with less than N unique "
                      "aligned fragments [default=%default]")
    parser.add_option("--isoform-fraction", type="float", 
                      default=0.10, metavar="X",
                      help="Filter chimeras with expression ratio "
                      " less than X (0.0-1.0) relative to the wild-type "
                      "5' transcript level [default=%default]")
    parser.add_option("--false-pos", dest="false_pos_file",
                      default=None, 
                      help="File containing known false positive "
                      "transcript pairs to subtract from output")
    options, args = parser.parse_args()
    index_dir = args[0]
    bam_file = args[1]
    input_file = args[2]
    output_file = args[3]
    return filter_chimeras(input_file, output_file, index_dir, bam_file,
                           unique_frags=options.unique_frags,
                           isoform_fraction=options.isoform_fraction,
                           false_pos_file=options.false_pos_file)

if __name__ == "__main__":
    main()