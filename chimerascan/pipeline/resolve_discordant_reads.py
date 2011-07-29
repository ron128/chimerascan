'''
Created on Jul 28, 2011

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

from chimerascan.lib.base import make_temp
from chimerascan.lib.chimera import Chimera
from chimerascan.pipeline.profile_insert_size import InsertSizeDistribution
from chimerascan.lib.batch_sort import batch_sort

QNAME_COL = 0
CHIMERA_NAME_COL = 1
SPANNING_FRAGS_COL = 2
NUM_UNAMBIGUOUS_FRAGS_COL = 3
NUM_UNIQUELY_ALIGNING_FRAGS_COL = 4
NUM_MISMATCHES_COL = 5
ISIZE_PROB_COL = 6
            
def choose_most_likely_alignments(dreads, chimeras):
    # criteria for selection the best alignment includes:
    # - total reads supporting chimera
    # - chimera with best spanning reads
    # - number of mismatches
    # - implied insert size of alignment
    pass

def merge_mmap_hists(hist1, hist2):
    newhist = []
    for count1,count2 in zip(hist1, hist2):
        newhist.append(count1 + count2)
    return newhist


def calc_isize_prob(isize, isize_dist):
    # find percentile of observing this insert size in the reads
    isize_per = isize_dist.percentile_at_isize(isize)
    # convert to a probability score (0.0-1.0)
    isize_prob = 1.0 - (2.0 * abs(50.0 - isize_per))/100.0    
    return isize_prob

def parse_read_stats(line_iter):
    for line in line_iter:
        fields = line.strip().split('\t')
        fields[2:6] = map(int, fields[2:6])
        fields[6] = float(fields[6])
        yield fields
        
def group_by_field(item_iter, colnum):
    mylist = []
    prev = None
    for fields in item_iter:
        # parse read stats information
        cur = fields[colnum]
        if prev != cur:
            if len(mylist) > 0:
                yield prev, mylist
                mylist = []
            prev = cur
        mylist.append(fields)
    if len(mylist) > 0:
        yield prev, mylist

def parse_sync_chimeras_read_stats(chimera_file, read_stats_file):
    # group reads by chimera name
    read_stats_iter = group_by_field(parse_read_stats(open(read_stats_file)), 1)
    iter_valid = True
    try:
        read_chimera_name, stats = read_stats_iter.next()
    except StopIteration:
        iter_valid = False
        stats = []
    # group chimeras by name    
    for c in Chimera.parse(open(chimera_file)):        
        while (iter_valid) and (c.name > read_chimera_name):
            try:
                read_chimera_name, stats = read_stats_iter.next()
            except StopIteration:
                iter_valid = False
                stats = []
        if c.name < read_chimera_name:
            yield c, []
        else:
            yield c, stats

def resolve_discordant_reads(input_file, output_file, isize_dist, min_isize_prob,
                             tmp_dir):
    #
    # parse chimeras and output reads to a file
    #
    logging.debug("Getting discordant read information")
    read_stats_file = os.path.join(tmp_dir, "read_stats.txt")
    f = open(read_stats_file, "w")
    for c in Chimera.parse(open(input_file)):
        # get number of unique alignment positions
        num_uniquely_aligning_frags = c.get_num_unique_positions()
        # get number of unambiguous reads
        num_unambiguous_frags = c.get_num_frags(maxnumhits=1)
        # number of spanning frags
        num_spanning_frags = c.get_num_spanning_frags()             
        # TODO: some statistics about spanning?
        for dpair in c.encomp_frags:
            # get putative insert size
            isize5p = c.tx_end_5p - dpair[0].pos
            isize3p = dpair[1].pos - c.tx_start_3p
            isize = isize5p + isize3p
            isize_prob = calc_isize_prob(isize, isize_dist)
            # output to file
            print >>f, '\t'.join(map(str, [dpair[0].qname, # read name 
                                           c.name, # chimera name
                                           num_spanning_frags,
                                           num_unambiguous_frags,
                                           num_uniquely_aligning_frags,
                                           -(dpair[0].mismatches + dpair[1].mismatches), # total mismatches
                                           isize_prob # insert size probability
                                           ]))
    f.close()
    #
    # now sort the read/chimera stats list
    #
    logging.debug("Sorting reads by read name")
    def sort_read_name(line):
        return line.strip().split('\t', 1)[QNAME_COL]
    sorted_read_stats_file = os.path.join(tmp_dir, "read_stats.rname_sorted.txt")
    batch_sort(input=read_stats_file,
               output=sorted_read_stats_file,
               key=sort_read_name,
               buffer_size=32000,
               tempdirs=[tmp_dir])
    #
    # parse reads by read name
    #
    logging.debug("Choosing best read groups")
    resolved_read_stats_file = os.path.join(tmp_dir, "read_stats.rname_sorted.resolved.txt")
    f = open(resolved_read_stats_file, "w")
    for rname,readstats in group_by_field(parse_read_stats(open(sorted_read_stats_file)), QNAME_COL):
        # build a dictionary of stats -> read/chimeras
        stats_dict = collections.defaultdict(lambda: [])
        for fields in readstats:
            # make a key to sort on
            key = tuple(fields[2:7])
            # value is chimera name
            stats_dict[key].append(fields[CHIMERA_NAME_COL])
        # sort based on stats
        sorted_stats_keys = sorted(stats_dict.keys(), reverse=True)
        best_key = sorted_stats_keys[0]
        # use only the best key
        chimera_names = stats_dict[best_key]
        # output read -> chimera relationships
        for chimera_name in chimera_names:
            print >>f, '\t'.join([rname, chimera_name] + map(str,best_key)) 
    f.close()
    #
    # re-sort by chimera name
    #
    logging.debug("Resorting reads by chimera name")
    def sort_reads_by_chimera_name(line):
        return line.strip().split('\t', 2)[CHIMERA_NAME_COL]
    sorted_resolved_read_stats_file = os.path.join(tmp_dir, "read_stats.chimera_name_sorted.resolved.txt")
    batch_sort(input=resolved_read_stats_file,
               output=sorted_resolved_read_stats_file,
               key=sort_reads_by_chimera_name,
               buffer_size=32000,
               tempdirs=[tmp_dir])
    logging.debug("Resorting chimeras by name")
    def sort_chimeras_by_name(line):
        return line.strip().split('\t', 7)[Chimera.NAME_FIELD]
    sorted_chimera_file = os.path.join(tmp_dir, "spanning_chimeras.name_sorted.txt")
    batch_sort(input=input_file,
               output=sorted_chimera_file,
               key=sort_chimeras_by_name,
               buffer_size=32000,
               tempdirs=[tmp_dir])
    #
    # parse and rebuild chimeras based on best reads
    # 
    logging.debug("Rewriting chimeras with lists of 'best' reads")
    f = open(output_file, "w")
    # need to sync chimeras with stats
    for c,stats in parse_sync_chimeras_read_stats(sorted_chimera_file, sorted_resolved_read_stats_file):
        # replace encompassing frags with reads in stats
        good_qnames = set(x[QNAME_COL] for x in stats
                          if x[ISIZE_PROB_COL] >= min_isize_prob)
        new_encomp_frags = [dpair for dpair in c.encomp_frags
                            if dpair[0].qname in good_qnames]
        c.encomp_frags = new_encomp_frags
        print >>f, '\t'.join(map(str, c.to_list()))
    f.close()
    # remove temporary files
    #os.remove(read_stats_file)
    #os.remove(sorted_read_stats_file)
    #os.remove(resolved_read_stats_file)
    #os.remove(sorted_resolved_read_stats_file)
    #os.remove(sorted_chimera_file)

def main():
    from optparse import OptionParser
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = OptionParser("usage: %prog [options] <in.txt> <out.txt> <isizedist.txt>")
    parser.add_option("--min-isize-prob", dest="min_isize_prob", 
                      type="float", default=0.01)
    options, args = parser.parse_args()
    input_file = args[0]
    output_file = args[1]
    isize_dist_file = args[2]
    # read insert size distribution
    isize_dist = InsertSizeDistribution()
    isize_dist.from_file(open(isize_dist_file, "r"))
    resolve_discordant_reads(input_file, output_file, isize_dist, 
                             options.min_isize_prob,
                             tmp_dir="/tmp")


if __name__ == '__main__':
    main()