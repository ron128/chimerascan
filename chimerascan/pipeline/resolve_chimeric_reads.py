'''
Created on Jul 1, 2011

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

from chimerascan.lib.chimera import Chimera
from chimerascan.pipeline.profile_insert_size import InsertSizeDistribution

class ChimeraStats(object):
    """
    stats used to assess quality of chimera
    """
    pass

def choose_most_likely_alignments(dreads, chimeras):
    # criteria for selection the best alignment includes:
    # - total reads supporting chimera
    # - chimera with best spanning reads
    # - number of mismatches
    # - implied insert size of alignment
    pass

def calc_isize_prob(dpair, c, isize_dist):
    # calculate insert size of read
    isize5p = c.partner5p.end - dpair[0].pos
    isize3p = dpair[1].aend - c.partner3p.start
    # find percentile of observing this insert size in the reads
    isize_per = isize_dist.percentile_at_isize(isize5p + isize3p)
    # convert to a probability score (0.0-1.0)
    isize_prob = 1.0 - (2.0 * abs(50.0 - isize_per))/100.0
    return isize_prob

def merge_mmap_hists(hist1, hist2):
    newhist = []
    for count1,count2 in zip(hist1, hist2):
        newhist.append(count1 + count2)
    return newhist

def resolve_chimeric_reads(input_file, output_file, isize_dist):
    # gather statistics on read alignments and chimeras
    # that will be used to associate reads with chimeras
    read_chimera_dict = collections.defaultdict(lambda: [])
    chimera_stats_dict = {}
    for c in Chimera.parse(open(input_file)):
        # combine multimap hists
        mmap_hist = merge_mmap_hists(c.partner5p.multimap_hist, 
                                     c.partner3p.multimap_hist)
        chimera_stats_dict[c.name] = (c.get_total_unique_reads(), 
                                      c.get_unique_spanning_reads(),
                                      mmap_hist)
        # get statistics for reads
        for dpair in c.encomp_read_pairs:
            qname = dpair[0].qname
            # find insert size probability of read
            isize_prob = calc_isize_prob(dpair, c, isize_dist)
            # find number of mismatches
            mismatches = dpair[0].mismatches + dpair[1].mismatches
            read_chimera_dict[qname].append((c.name, isize_prob, mismatches))
    # now process one read at a time, looking at all its alignments and
    # choose the most likely set of alignments
    encomp_qnames_dict = collections.defaultdict(lambda: [])
    for qname, read_stats_list in read_chimera_dict.iteritems():
        stats_chimera_dict = collections.defaultdict(lambda: [])
        for chimera_name, isize_prob, mismatches in read_stats_list:
            # get chimera stats
            unique_reads, spanning_reads, mmap_hist = chimera_stats_dict[chimera_name]
            # make a key to sort on
            key = (unique_reads, spanning_reads) + tuple(mmap_hist) + (isize_prob, -mismatches)
            stats_chimera_dict[key].append(chimera_name)
        # sort keys (reverse)
        sorted_stats_keys = sorted(stats_chimera_dict.keys(), reverse=True)
        # use only the best key
        chimera_names = stats_chimera_dict[sorted_stats_keys[0]]
        for chimera_name in chimera_names:
            encomp_qnames_dict[chimera_name].add(qname)
    # now edit the chimeras using the modified qnames
    f = open(output_file, "w")
    for c in Chimera.parse(open(input_file)):
        qnames = encomp_qnames_dict[c.name]
        filtered_encomp_pairs = []
        for pair in c.encomp_read_pairs:
            if pair[0].qname not in qnames:
                continue
            filtered_encomp_pairs.append(pair)
        # update encompassing reads
        c.encomp_read_pairs = filtered_encomp_pairs
        
    f.close()


def main():
    from optparse import OptionParser
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = OptionParser("usage: %prog [options] <in.txt> <out.txt> <isizedist.txt>")
    parser.add_option("--total-reads", type="int", default=2,
                      dest="total_reads_threshold", metavar="N",
                      help="Filter chimeras with less than N total "
                      "unique supporting reads")
    options, args = parser.parse_args()
    input_file = args[0]
    output_file = args[1]
    isize_dist_file = args[2]
    # read insert size distribution
    isize_dist = InsertSizeDistribution()
    isize_dist.from_file(open(isize_dist_file, "r"))
    resolve_chimeric_reads(input_file, output_file, isize_dist)


if __name__ == '__main__':
    main()