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
import numpy as np

# local imports
from chimerascan import pysam
from chimerascan.lib.alignment_parser import parse_sr_sam_file
from nominate_chimeras import Chimera

def read_chimera_mapping_file(filename):
    chimera_refs = collections.defaultdict(lambda: [])
    for line in open(filename):
        fields = line.strip().split('\t')
        chimera_id = fields[0]
        bedpe_fields = fields[1:]        
        chimera_refs[chimera_id].append(bedpe_fields)
    return dict(chimera_refs)

def get_mismatch_positions(md):
    x = 0
    pos = []
    for y in xrange(len(md)):
        if md[y].isalpha():
            offset = int(md[x:y])
            pos.append(offset)
            x = y + 1
    return pos

def filter_reads_by_anchor(reads, junc_positions, anchors,
                           anchor_min, anchor_max, max_anchor_mismatches):
    for read in reads:
        if read.is_unmapped:
            continue
        # filter out reads with small anchor sequence
        aligned_length = read.aend - read.pos  
        junc_pos = junc_positions[read.rname]
        # check if the read spans the junction
        read_ok = True
        if read.aend <= junc_pos or read.pos >= junc_pos:
            # does not span
            read_ok = False
        else:
            left_anchor_bp = junc_pos - read.pos
            right_anchor_bp = read.aend - junc_pos
            # keep track of minimum number of bp spanning the chimera
            anchor = min(left_anchor_bp, right_anchor_bp)
            anchors[read.rname][anchor - 1] += 1
            if anchor < anchor_min:
                # not enough anchor
                read_ok = False
            elif anchor < anchor_max:
                # find anchor interval
                if left_anchor_bp < anchor_max:
                    anchor_interval = (0, left_anchor_bp)
                else:
                    anchor_interval = (aligned_length - right_anchor_bp, aligned_length)      
                # get mismatch positions
                mmpos = get_mismatch_positions(read.opt('MD'))
                # see if mismatches lie in anchor interval
                anchor_mm = [pos for pos in mmpos
                             if pos >= anchor_interval[0] and pos < anchor_interval[1]]
                if len(anchor_mm) > max_anchor_mismatches:
                    # mismatches within anchor position
                    read_ok = False
        if read_ok:
            yield read

def join_spanning_reads(bam_file, junc_map,
                        original_read_length,
                        anchor_min, anchor_max,
                        max_anchor_mismatches):
    # map reference names to numeric ids
    bamfh = pysam.Samfile(bam_file, "rb")
    rname_tid_map = dict((rname,i) for i,rname in enumerate(bamfh.references))
    anchors = {}
    junc_positions = {}
    max_anchor = int((original_read_length + 1)/2)
    for chimera_name in junc_map.iterkeys():
        left_junc_length = int(chimera_name.split(":")[1])
        chimera_id = rname_tid_map[chimera_name]
        junc_positions[chimera_id] = left_junc_length
        anchors[chimera_id] = np.zeros(max_anchor, dtype=np.int)
    # count read alignments to chimera junctions
    chimera_counts = collections.defaultdict(lambda: 0)
    chimera_reads = collections.defaultdict(lambda: [])
    num_multimaps = 0    
    num_filtered = 0
    num_reads = 0
    for reads in parse_sr_sam_file(bamfh):
        num_reads += 1
        if len(reads) > 1:
            num_multimaps += 1
        filtered_reads = \
            list(filter_reads_by_anchor(reads, junc_positions, anchors,
                                        anchor_min, anchor_max, max_anchor_mismatches))
        if len(filtered_reads) == 0:
            # no reads passed filter
            num_filtered += 1
            continue
        for read in filtered_reads:
            chimera_counts[read.rname] += 1
            chimera_reads[read.rname].append((read.qname,read.seq))
    logging.info("Reads: %d" % (num_reads))
    logging.info("Multimapping: %d" % (num_multimaps))
    logging.info("Failed anchor filter: %d" % (num_filtered))
    bamfh.close()
    # assign junction coverage to chimera candidates
    for chimera_name,bedpe_records in junc_map.iteritems():
        chimera_id = rname_tid_map[chimera_name]        
        cov = chimera_counts[chimera_id]
        read_tuples = chimera_reads[chimera_id]
        spanning_qnames = set(read_tuple[0] for read_tuple in read_tuples)
        if len(read_tuples) == 0:
            read_tuples = [("None", "None")]
        for fields in bedpe_records:
            # intersect encompassing with spanning reads to see overlap
            encompassing_qnames = fields[Chimera.QNAME_COL].split(Chimera.SEQ_FIELD_DELIM)
            union_cov = len(spanning_qnames.union(encompassing_qnames))
            intersect_cov = len(spanning_qnames.intersection(encompassing_qnames))
            #both_cov = sum(1 for read_tuple in read_tuples if read_tuple[0] in set(encompassing_qnames))
            fields += [cov, intersect_cov, union_cov, ','.join(map(str,anchors[chimera_id])),
                       Chimera.SEQ_FIELD_DELIM.join([read_tuple[1] for read_tuple in read_tuples]),
                       Chimera.SEQ_FIELD_DELIM.join([read_tuple[0] for read_tuple in read_tuples])]
            yield fields

def merge_spanning_alignments(bam_file, junc_map_file, output_file,
                              read_length, anchor_min, anchor_max,
                              anchor_mismatches):
    junc_map = read_chimera_mapping_file(junc_map_file)
    f = open(output_file, "w")
    for bedpe_fields in join_spanning_reads(bam_file, junc_map,
                                            read_length,
                                            anchor_min,
                                            anchor_max, 
                                            anchor_mismatches):
        print >>f, '\t'.join(map(str, bedpe_fields))
    f.close()

def main():
    from optparse import OptionParser
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")    
    parser = OptionParser("usage: %prog [options] <in.bam> <junc_map> <out.txt>")    
    parser.add_option("--rlen", type="int", dest="read_length")
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



## test uniformity of reads spanning junction
#observed_arr = np.array(np.round(positions[chimera_id],0), dtype=np.int)
#expected_arr = np.zeros(observed_arr.shape[0], dtype=np.int)
## find integer and fraction parts of expected        
#expected_int = int(observed_arr.mean())
#expected_remainder = observed_arr.sum() % observed_arr.shape[0]
## randomly assign the remaining counts
#remainder_positions = np.random.randint(0, observed_arr.shape[0],size=expected_remainder)
#expected_arr[:] = expected_int
#for pos in remainder_positions:
#    expected_arr[pos] += 1
## perform chi-squared test for coverage uniformity
#csq, pval = chisquare(observed_arr, expected_arr)


#    num_multiple_partitions = 0
#    num_splits = 0
#    num_filtered = 0
#    for alignments in parse_segmented_sr_sam_file(bamfh):
#        filtered_reads = []
#        if len(alignments) > 1:
#            num_multiple_partitions += 1
#        for partition in alignments:
#            if len(partition) > 1:
#                num_splits += 1
#                # TODO: skip reads that split into multiple partitions
#                # since the junctions should be contiguous
#                continue
#            for split_reads in partition:
#                func = filter_reads_by_anchor(split_reads, junc_positions, anchors,
#                                              anchor_min, anchor_max, max_anchor_mismatches)
#                filtered_reads.extend(func)
#        if len(filtered_reads) == 0:
#            # no reads passed filter
#            num_filtered += 1
#            continue
#        for read in filtered_reads:
#            chimera_counts[read.rname] += 1
#            chimera_reads[read.rname].append((read.qname,read.seq))