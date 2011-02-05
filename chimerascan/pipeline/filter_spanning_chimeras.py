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
import collections
import logging
import os

# local lib imports
from chimerascan.lib import config
from chimerascan.lib.gene_to_genome2 import build_gene_to_genome_map, gene_to_genome_pos
# local imports
from nominate_chimeras import Chimera
        
class SpanningChimera(Chimera):
    FIRST_COL = Chimera.LAST_COL + 1
    LAST_COL = FIRST_COL + 6 
    
    def __init__(self):
        Chimera.__init__(self)
        self.spanning_reads = 0
        self.encomp_and_spanning = 0
        self.total_reads = 0
        self.junction_hist = None
        self.junction_kl_div = 0.0
        self.spanning_seq = None
        self.spanning_ids = None
        
    def from_list(self, fields):
        FIRST_COL = Chimera.LAST_COL + 1
        # get the chimera fields
        Chimera.from_list(self, fields)
        self.spanning_reads = int(fields[FIRST_COL])
        self.encomp_and_spanning = int(fields[FIRST_COL+1])
        self.total_reads = int(fields[FIRST_COL+2])
        self.junction_hist = map(int, fields[FIRST_COL+3].split(','))
        self.junction_kl_div = float(fields[FIRST_COL+4])
        self.spanning_seq = fields[FIRST_COL+5]
        self.spanning_ids = fields[FIRST_COL+6]

    def to_list(self):
        fields = Chimera.to_list(self)
        fields.extend([self.spanning_reads, 
                       self.encomp_and_spanning, 
                       self.total_reads, 
                       ','.join(map(str, self.junction_hist)), 
                       self.junction_kl_div,
                       self.spanning_seq, 
                       self.spanning_ids])
        return fields

    @staticmethod
    def parse(line_iter):
        for line in line_iter:
            fields = line.strip().split('\t')
            c = SpanningChimera()
            c.from_list(fields)
            yield c

def build_junc_coverage_map(chimeras, ggmap):
    junc_cov_map = collections.defaultdict(lambda: [None, None, None])
    num_chimeras = 0
    for c in chimeras:
        num_chimeras += 1
        # convert to genomic coords
        # subtract one since 5' junc position is an open interval
        coord5p = gene_to_genome_pos(c.mate5p.tx_name, c.mate5p.end - 1, ggmap)
        coord3p = gene_to_genome_pos(c.mate3p.tx_name, c.mate3p.start, ggmap)
        # keep track of maximum coverage isoform
        pairkey = (coord5p, coord3p)
        paircov = (c.encomp_and_spanning, c.weighted_cov, c.total_reads)
        data = junc_cov_map[pairkey]
        if (data[0] is None) or (cmp(paircov, data[0]) > 0):
            # store encomp/spanning, weighted coverage, and total reads
            data[0] = paircov
            data[1] = c.mate5p.tx_name
            data[2] = c.mate3p.tx_name
    logging.debug("Parsed %d chimeras" % (num_chimeras))
    kept_isoforms = set(tuple(v[1:3]) for v in junc_cov_map.itervalues())
    #del junc_cov_map
    logging.debug("Kept %d highest coverage isoforms" % (len(kept_isoforms)))
    return kept_isoforms

def choose_highest_coverage_chimeras(input_file, ggmap):
    '''
    choose the highest coverage isoform pair using spanning reads,
    encompassing reads, and total reads as a measure.  ties will be
    broken by choosing a single gene pair arbitrarily 
    '''
    # break name into 5'/3' genes linked in a dictionary
    logging.debug("Building junction isoform coverage map")
    kept_isoforms_set = build_junc_coverage_map(SpanningChimera.parse(open(input_file)), ggmap)
    # write results
    logging.debug("Returning highest coverage chimeras")
    for c in SpanningChimera.parse(open(input_file)):
        pairkey = (c.mate5p.tx_name, c.mate3p.tx_name)
        if pairkey in kept_isoforms_set:
            yield c
    del kept_isoforms_set

def filter_spanning_chimeras(input_file, output_file, gene_file):
    logging.debug("Building gene/genome index")
    ggmap = build_gene_to_genome_map(open(gene_file))
    logging.debug("Choosing highest coverage chimeras")
    fh = open(output_file, "w")
    for c in choose_highest_coverage_chimeras(input_file, ggmap):
        print >>fh, '\t'.join(['\t'.join(map(str,c.to_list()))])
    fh.close()

def main():
    from optparse import OptionParser
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = OptionParser("usage: %prog [options] <sortedchimeras.bedpe> <chimeras.txt>")
    parser.add_option("--index", dest="index_dir",
                      help="Path to chimerascan index directory")
    options, args = parser.parse_args()
    gene_file = os.path.join(options.index_dir, config.GENE_FEATURE_FILE)
    input_file = args[0]
    output_file = args[1]
    filter_spanning_chimeras(input_file, output_file, gene_file)

if __name__ == "__main__":
    main()