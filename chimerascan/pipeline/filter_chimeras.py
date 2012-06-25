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
import sys
import os
import logging
import collections
import argparse

import pysam

from chimerascan.lib.transcriptome import build_transcript_genome_map, \
    transcript_to_genome_pos
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

def filter_chimeras(input_file, output_file,
                    filter_num_frags,
                    filter_allele_fraction,
                    mask_biotypes,
                    mask_rnames):
    logging.debug("Filtering chimeras")
    logging.debug("\tfragments: %f" % (filter_num_frags))
    logging.debug("\tallele fraction: %f" % (filter_allele_fraction))
    logging.debug("\tmask biotypes: %s" % (','.join(sorted(mask_biotypes))))
    logging.debug("\tmask references: %s" % (','.join(sorted(mask_rnames))))
    # filter chimeras
    num_chimeras = 0
    num_kept_chimeras = 0    
    f = open(output_file, "w")   
    for c in Chimera.parse(open(input_file)):
        num_chimeras += 1
        # number of fragments
        if c.num_frags < filter_num_frags:
            continue
        # allele fraction
        allele_fraction_5p = float(c.num_frags) / (c.num_discordant_frags_5p + c.num_concordant_frags_5p)
        allele_fraction_3p = float(c.num_frags) / (c.num_discordant_frags_3p + c.num_concordant_frags_3p)
        allele_fraction = min(allele_fraction_5p, allele_fraction_3p)
        if allele_fraction < filter_allele_fraction:
            continue
        # masked biotypes and references
        if len(mask_biotypes.intersection(c.biotypes_5p)) > 0:
            continue
        if len(mask_biotypes.intersection(c.biotypes_3p)) > 0:
            continue
        if c.rname5p in mask_rnames:
            continue
        if c.rname3p in mask_rnames:
            continue
        print >>f, str(c)
        num_kept_chimeras += 1
    f.close()
    logging.debug("Total chimeras: %d" % num_chimeras)
    logging.debug("Kept chimeras: %d" % num_kept_chimeras)
    return config.JOB_SUCCESS

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("--num-frags", dest="num_frags", 
                        type=float, default=config.DEFAULT_FILTER_FRAGS)
    parser.add_argument("--allele-fraction", type=float, 
                        default=config.DEFAULT_FILTER_ALLELE_FRACTION, 
                        dest="allele_fraction", metavar="X",
                        help="Filter chimeras with expression less than "
                        "the specified fraction of the total expression "
                        "level [default=%(default)s")
    parser.add_argument("--mask-biotypes-file", dest="mask_biotypes_file", default=None) 
    parser.add_argument("--mask-rnames-file", dest="mask_rnames_file", default=None)
    parser.add_argument("input_file")
    parser.add_argument("output_file")
    args = parser.parse_args()
    mask_biotypes = set()
    if args.mask_biotypes_file is not None:
        mask_biotypes.update([line.strip() for line in open(args.mask_biotypes_file)])
    mask_rnames = set()
    if args.mask_rnames_file is not None:
        mask_rnames.update([line.strip() for line in open(args.mask_rnames_file)])
    return filter_chimeras(args.input_file, args.output_file,
                           filter_num_frags=args.num_frags,
                           filter_allele_fraction=args.allele_fraction,
                           mask_biotypes=mask_biotypes,
                           mask_rnames=mask_rnames)

if __name__ == "__main__":
    sys.exit(main())