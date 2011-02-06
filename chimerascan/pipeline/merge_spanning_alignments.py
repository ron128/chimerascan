'''
Created on Nov 7, 2010

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
import operator

# local imports
from chimerascan import pysam
from chimerascan.lib.alignment_parser import parse_sr_sam_file
from chimerascan.lib.base import parse_string_none, select_best_mismatch_strata
from chimerascan.lib.stats import kl_divergence
from nominate_chimeras import Chimera

SpanningRead = collections.namedtuple("SpanningRead", 
                                      ["qname", "mate", "seq", "pos", 
                                       "aend", "mappings", "is_reverse"])

class SpanningChimera(Chimera):
    FIRST_COL = Chimera.LAST_COL + 1
    LAST_COL = FIRST_COL + 5     
    SPAN_READ_DELIM = ';'
    SPAN_FIELD_DELIM = '|'
    
    def __init__(self):
        Chimera.__init__(self)
        self.junc_name = None
        self.junc_pos = 0
        self.junc_homology_5p = 0
        self.junc_homology_3p = 0
        self.num_spanning_reads = 0
        self.encomp_and_spanning = 0
        self.encomp_or_spanning = 0
        self.spanning_reads = []
        
    def from_list(self, fields):
        FIRST_COL = Chimera.LAST_COL + 1
        # get the chimera fields
        Chimera.from_list(self, fields)
        self.junc_name = fields[FIRST_COL]
        self.junc_pos = int(fields[FIRST_COL+1])
        self.junc_homology_5p = int(fields[FIRST_COL+2])
        self.junc_homology_3p = int(fields[FIRST_COL+3])        
        self.num_spanning_reads = int(fields[FIRST_COL+4])
        self.encomp_and_spanning = int(fields[FIRST_COL+5])
        self.encomp_or_spanning = int(fields[FIRST_COL+6])        
        spanning_reads_field = parse_string_none(fields[FIRST_COL+7]) 
        self.spanning_reads = []
        if spanning_reads_field is not None:
            for read_fields in spanning_reads_field.split(self.SPAN_READ_DELIM):
                fields = read_fields.split(self.SPAN_FIELD_DELIM)
                read = SpanningRead(qname=fields[0], 
                                    mate=int(fields[1]),
                                    seq=fields[2],
                                    pos=int(fields[3]), 
                                    aend=int(fields[4]),
                                    mappings=int(fields[5]),
                                    is_reverse=int(fields[6]))
                self.spanning_reads.append(read)

    def to_list(self):
        fields = Chimera.to_list(self)[:Chimera.LAST_COL+1]
        fields.extend([self.junc_name,
                       self.junc_pos,
                       self.junc_homology_5p,
                       self.junc_homology_3p,
                       self.num_spanning_reads,
                       self.encomp_and_spanning,
                       self.encomp_or_spanning])   
        read_fields = [self.SPAN_FIELD_DELIM.join(map(str, read)) 
                       for read in self.spanning_reads]
        if len(read_fields) == 0:
            fields.append("None")
        else:
            fields.append(self.SPAN_READ_DELIM.join(read_fields))
        return fields

    @staticmethod
    def parse(line_iter):
        for line in line_iter:
            fields = line.strip().split('\t')
            c = SpanningChimera()
            c.from_list(fields)
            yield c

class SpanningReads(object):
    def __init__(self):
        self.junc_tid = None
        self.junc_pos = None
        self.junc_homology_5p = 0
        self.junc_homology_5p = 0
        self.reads = []

def read_junc_mapping_file(junc_map_fh, rname_tid_map):
    spanning_data = collections.defaultdict(lambda: SpanningReads())
    for line in junc_map_fh:
        fields = line.strip().split('\t')
        junc_tid = rname_tid_map[fields[0]]        
        s = spanning_data[junc_tid]
        s.junc_tid = junc_tid
        s.junc_pos = int(fields[1])
        s.junc_homology_5p = int(fields[2])
        s.junc_homology_3p = int(fields[3])
    return spanning_data

def get_mismatch_positions(md):
    x = 0
    pos = []
    for y in xrange(len(md)):
        if md[y].isalpha():
            offset = int(md[x:y])
            pos.append(offset)
            x = y + 1
    return pos

def filter_anchor_position(read,
                           spanning_data,
                           anchor_min,
                           anchor_max,
                           max_anchor_mismatches):
    junc_pos = spanning_data.junc_pos
    # check if the read spans the junction
    passes_filter = read.pos < junc_pos < read.aend
    if passes_filter:
        # read spans junction, but might violate anchor constraints        
        left_anchor_bp = junc_pos - read.pos
        right_anchor_bp = read.aend - junc_pos
        # check 5' homology
        if left_anchor_bp <= (spanning_data.junc_homology_5p + anchor_min):
            #logging.debug("Failed 5' homology filter left anchor=%d homology 5p=%d" %
            #              (left_anchor_bp, spanning_data.junc_homology_5p))
            passes_filter = False
        # check 3' homology
        if right_anchor_bp <= (spanning_data.junc_homology_3p + anchor_min):
            #logging.debug("Failed 3' homology filter right anchor=%d homology 3p=%d" %
            #              (right_anchor_bp, spanning_data.junc_homology_3p))
            passes_filter = False
        # keep track of minimum number of bp spanning the chimera
        anchor = min(left_anchor_bp, right_anchor_bp)
        if anchor < anchor_max:
            # count mismatches in anchor region
            # find anchor interval
            if left_anchor_bp < anchor_max:
                anchor_interval = (0, left_anchor_bp)
            else:
                aligned_length = read.aend - read.pos
                anchor_interval = (aligned_length - right_anchor_bp, aligned_length)      
            # get mismatch positions
            mmpos = get_mismatch_positions(read.opt('MD'))
            # see if any mismatches lie in anchor interval
            anchor_mm = [pos for pos in mmpos
                         if pos >= anchor_interval[0] and pos < anchor_interval[1]]
            if len(anchor_mm) > max_anchor_mismatches:
                # mismatches within anchor position
                passes_filter = False
    return passes_filter

def process_spanning_reads(reads, 
                           spanning_data_dict,
                           anchor_min=0,
                           anchor_max=0,
                           max_anchor_mismatches=0):
    passed_filter = False
    num_mapped_reads = sum(1 for r in reads if not r.is_unmapped)
    for read in reads:
        if read.is_unmapped:
            continue
        # get spanning read information
        spandata = spanning_data_dict[read.rname]
        # check if the read spans the junction and
        # passes anchor constraints
        if filter_anchor_position(read, spandata, anchor_min,
                                  anchor_max, max_anchor_mismatches):        
            # add read information to dict
            spandata.reads.append(SpanningRead(qname=read.qname,
                                               mate=int(read.is_read2), 
                                               seq=read.seq, 
                                               pos=read.pos, 
                                               aend=read.aend,
                                               mappings=num_mapped_reads, 
                                               is_reverse=int(read.is_reverse)))
            passed_filter = True
    return passed_filter


def parse_spanning_bam(bam_file,
                       junc_mapping_file,
                       anchor_min, anchor_max,
                       max_anchor_mismatches):
    # map reference names to numeric ids
    bamfh = pysam.Samfile(bam_file, "rb")
    rname_tid_map = dict((rname,i) for i,rname in enumerate(bamfh.references))
    spanning_data_dict = read_junc_mapping_file(open(junc_mapping_file),
                                                rname_tid_map)
    # count read alignments to chimera junctions
    num_multimaps = 0    
    num_filtered = 0
    num_reads = 0    
    for reads in parse_sr_sam_file(bamfh):
        num_reads += 1
        if len(reads) > 1:
            num_multimaps += 1
        passed_filter = process_spanning_reads(reads, spanning_data_dict, 
                                               anchor_min, anchor_max, 
                                               max_anchor_mismatches)
        if not passed_filter:
            # no reads passed filter
            num_filtered += 1
    logging.info("Reads: %d" % (num_reads))
    logging.info("Multimapping: %d" % (num_multimaps))
    logging.info("Failed anchor filter: %d" % (num_filtered))
    bamfh.close()
    return rname_tid_map, spanning_data_dict

def make_spanning_chimeras(spanning_data_dict, junc_mapping_file, rname_tid_map):    
    # output chimera candidates
    for line in open(junc_mapping_file):
        fields = line.strip().split('\t')
        # create spanning chimera object with all data 
        c = SpanningChimera()
        # first field in the junction map is the junction id
        junc_name = fields[0]        
        junc_tid = rname_tid_map[junc_name]        
        # fill in junction map fields
        c.junc_name = junc_name
        c.junc_pos = int(fields[1])
        c.junc_homology_5p = int(fields[2])
        c.junc_homology_3p = int(fields[3])
        # initialize with chimera fields
        Chimera.from_list(c, fields[4:])
        # spanning data
        spandata = spanning_data_dict[junc_tid]
        c.spanning_reads = spandata.reads
        # compute statistics
        spanning_qnames = set(r.qname for r in spandata.reads)        
        c.num_spanning_reads = len(spanning_qnames)
        c.encomp_and_spanning = len(spanning_qnames.intersection(c.qnames))
        c.encomp_or_spanning = len(spanning_qnames.union(c.qnames))
        yield c

def merge_spanning_alignments(bam_file, junc_map_file, output_file,
                              anchor_min, anchor_max,
                              anchor_mismatches):
    f = open(output_file, "w")
    rname_tid_map, spanning_data_dict = \
        parse_spanning_bam(bam_file, junc_map_file, anchor_min, anchor_max,
                           anchor_mismatches)
    for c in make_spanning_chimeras(spanning_data_dict, junc_map_file, 
                                    rname_tid_map):    
        print >>f, '\t'.join(map(str, c.to_list()))
    f.close()

def main():
    from optparse import OptionParser
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")    
    parser = OptionParser("usage: %prog [options] <in.bam> <junc_map> <out.txt>")    
    parser.add_option("--anchor-min", type="int", dest="anchor_min", default=0)
    parser.add_option("--anchor-max", type="int", dest="anchor_max", default=0)
    parser.add_option("--anchor-mismatches", type="int", dest="anchor_mismatches", default=0)
    options, args = parser.parse_args()
    bam_file = args[0]
    junc_map_file = args[1]
    output_file = args[2]
    merge_spanning_alignments(bam_file, junc_map_file, output_file,
                              options.read_length, 
                              options.anchor_min, 
                              options.anchor_max,
                              options.anchor_mismatches)

if __name__ == '__main__': main()
