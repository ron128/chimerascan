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

from chimerascan.lib.chimera import Chimera
from chimerascan.pipeline.profile_insert_size import InsertSizeDistribution
from chimerascan.lib.batch_sort import batch_sort

QNAME_COL = 0
CHIMERA_NAME_COL = 5
SCORE_FIELDS = (6,7,8,9,10)

class ChimeraStats(object):
    __slots__ = ('qname', 'tid5p', 'pos5p', 'tid3p', 'pos3p', 
                 'chimera_name', 'num_spanning_frags', 'num_unambiguous_frags',
                 'num_uniquely_aligning_frags', 'neg_mismatches',
                 'isize_prob')

    @property
    def score_tuple(self):
        return (self.num_spanning_frags,
                self.num_unambiguous_frags,
                self.num_uniquely_aligning_frags,
                self.neg_mismatches,
                self.isize_prob)

    def to_list(self):
        return [self.qname,
                self.tid5p, self.pos5p,
                self.tid3p, self.pos3p, 
                self.chimera_name,
                self.num_spanning_frags,
                self.num_unambiguous_frags,
                self.num_uniquely_aligning_frags,
                self.neg_mismatches,
                self.isize_prob]

    @staticmethod
    def from_list(fields):
        s = ChimeraStats()
        s.qname = fields[0]
        s.tid5p = int(fields[1])
        s.pos5p = int(fields[2])
        s.tid3p = int(fields[3])
        s.pos3p = int(fields[4])
        s.chimera_name = fields[5]
        s.num_spanning_frags = int(fields[6])
        s.num_unambiguous_frags = int(fields[7])
        s.num_uniquely_aligning_frags = int(fields[8])
        s.neg_mismatches = int(fields[9])
        s.isize_prob = float(fields[10])
        return s

    @staticmethod
    def parse(line_iter):
        for line in line_iter:
            fields = line.strip().split('\t')
            yield ChimeraStats.from_list(fields)

def calc_isize_prob(isize, isize_dist):
    # find percentile of observing this insert size in the reads
    isize_per = isize_dist.percentile_at_isize(isize)
    # convert to a probability score (0.0-1.0)
    isize_prob = 1.0 - (2.0 * abs(50.0 - isize_per))/100.0    
    return isize_prob

def group_by_attr(item_iter, attr):
    mylist = []
    prev = None
    for itm in item_iter:
        cur = getattr(itm, attr)
        if prev != cur:
            if len(mylist) > 0:
                yield prev, mylist
                mylist = []
            prev = cur
        mylist.append(itm)
    if len(mylist) > 0:
        yield prev, mylist

#def group_by_field(item_iter, colnum):
#    mylist = []
#    prev = None
#    for fields in item_iter:
#        # parse read stats information
#        cur = fields[colnum]
#        if prev != cur:
#            if len(mylist) > 0:
#                yield prev, mylist
#                mylist = []
#            prev = cur
#        mylist.append(fields)
#    if len(mylist) > 0:
#        yield prev, mylist

def parse_sync_chimeras_read_stats(chimera_file, read_stats_file):
    # group reads by chimera name
    read_stats_iter = group_by_attr(ChimeraStats.parse(open(read_stats_file)), 
                                    'chimera_name')
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
        for dpair in c.encomp_frags:
            # get putative insert size
            isize5p = c.tx_end_5p - dpair[0].pos
            isize3p = dpair[1].pos - c.tx_start_3p
            isize = isize5p + isize3p
            isize_prob = calc_isize_prob(isize, isize_dist)
            # make ChimeraStats object
            s = ChimeraStats()
            s.qname = dpair[0].qname
            s.tid5p = dpair[0].tid
            s.pos5p = dpair[0].pos
            s.tid3p = dpair[1].tid
            s.pos3p = dpair[1].pos
            s.chimera_name = c.name
            s.num_spanning_frags = num_spanning_frags
            s.num_unambiguous_frags = num_unambiguous_frags
            s.num_uniquely_aligning_frags = num_uniquely_aligning_frags
            s.neg_mismatches = -(dpair[0].mismatches + dpair[1].mismatches)
            s.isize_prob = isize_prob
            # output to file
            print >>f, '\t'.join(map(str, s.to_list()))
    f.close()
    #
    # now sort the read/chimera stats list
    #
    logging.debug("Sorting reads by read name")
    def sort_read_name(line):
        return line.strip().split('\t', QNAME_COL+1)[QNAME_COL]
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
    for rname,readstats in group_by_attr(ChimeraStats.parse(open(sorted_read_stats_file)), 
                                         'qname'):
        # build a dictionary of stats -> read/chimeras
        stats_dict = collections.defaultdict(lambda: [])
        for s in readstats:
            # add key/value pairs
            stats_dict[s.score_tuple].append(s)
        # sort based on stats
        sorted_stats_keys = sorted(stats_dict.keys(), reverse=True)
        # use only the best key
        for s in stats_dict[sorted_stats_keys[0]]:
            # output read -> chimera relationships
            print >>f, '\t'.join(map(str, s.to_list()))
    f.close()
    #
    # re-sort by chimera name
    #
    logging.debug("Resorting reads by chimera name")
    def sort_reads_by_chimera_name(line):
        return line.strip().split('\t',CHIMERA_NAME_COL+1)[CHIMERA_NAME_COL]
    sorted_resolved_read_stats_file = os.path.join(tmp_dir, "read_stats.chimera_name_sorted.resolved.txt")
    batch_sort(input=resolved_read_stats_file,
               output=sorted_resolved_read_stats_file,
               key=sort_reads_by_chimera_name,
               buffer_size=32000,
               tempdirs=[tmp_dir])
    logging.debug("Resorting chimeras by name")
    def sort_chimeras_by_name(line):
        return line.strip().split('\t',Chimera.NAME_FIELD+1)[Chimera.NAME_FIELD]
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
        # parse and make lookup set of the resolved alignments
        good_alignments = set()
        for s in stats:
            if s.isize_prob < min_isize_prob:
                continue
            good_alignments.add((s.qname, s.tid5p, s.pos5p, s.tid3p, s.pos3p))
        # replace encompassing frags with resolved alignments
        new_encomp_frags = []
        for dpair in c.encomp_frags:
            # get alignment tuple
            aln = (dpair[0].qname, dpair[0].tid, dpair[0].pos, dpair[1].tid, dpair[1].pos)
            if aln in good_alignments:
                new_encomp_frags.append(dpair)
        c.encomp_frags = new_encomp_frags
        c.score = c.get_num_frags()
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
    isize_dist = InsertSizeDistribution.from_file(open(isize_dist_file))
    resolve_discordant_reads(input_file, output_file, isize_dist, 
                             options.min_isize_prob,
                             tmp_dir=".")

if __name__ == '__main__':
    main()