'''
Created on Nov 7, 2010

@author: mkiyer
'''
import logging
import collections
import numpy as np

import pysam

# local imports

# constants
JUNC_MAP_QNAME_COLUMN = 14
SEQ_FIELD_DELIM = ';'

def read_chimera_mapping_file(filename):
    chimera_refs = collections.defaultdict(lambda: [])
    for line in open(filename):
        fields = line.strip().split('\t')
        chimera_id = fields[0]
        bedpe_fields = fields[1:]        
        chimera_refs[chimera_id].append(bedpe_fields)
    return dict(chimera_refs)

def parse_sr_sam_file(bamfh):
    reads = []
    # reads must be binned by qname, mate, hit, and segment
    # so initialize to mate 0, hit 0, segment 0
    num_reads = 0
    prev_qname = None
    for read in bamfh:
        # get read attributes
        qname = read.qname
        mate = 0 if read.is_read1 else 1
        # get hit/segment/mapping tags
        num_split_partitions = read.opt('NH')
        partition_ind = read.opt('XH')
        num_splits = read.opt('XN')
        split_ind = read.opt('XI')
        num_mappings = read.opt('IH')
        mapping_ind = read.opt('HI')        
        # if query name changes we have completely finished
        # the fragment and can reset the read data
        if num_reads > 0 and qname != prev_qname:
            yield reads
            # reset state variables
            reads = []
            num_reads = 0
        prev_qname = qname
        # initialize hits
        if len(reads) == 0:
            reads.extend([list() for x in xrange(num_split_partitions)])
        # initialize hit segments
        if len(reads[partition_ind]) == 0:
            reads[partition_ind].extend([list() for x in xrange(num_splits)])
        split_reads = reads[partition_ind][split_ind]
        # initialize segment mappings
        if len(split_reads) == 0:
            split_reads.extend([None for x in xrange(num_mappings)])
        # add segment to hit/mate/read
        split_reads[mapping_ind] = read
        num_reads += 1
    if num_reads > 0:
        yield reads

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
    num_multiple_partitions = 0
    num_splits = 0
    num_filtered = 0
    for alignments in parse_sr_sam_file(bamfh):
        filtered_reads = []
        if len(alignments) > 1:
            num_multiple_partitions += 1
        for partition in alignments:
            if len(partition) > 1:
                num_splits += 1
            for split_reads in partition:
                func = filter_reads_by_anchor(split_reads, junc_positions, anchors,
                                              anchor_min, anchor_max, max_anchor_mismatches)
                filtered_reads.extend(func)
        if len(filtered_reads) == 0:
            # no reads passed filter
            num_filtered += 1
            continue
        for read in filtered_reads:
            chimera_counts[read.rname] += 1
            chimera_reads[read.rname].append((read.qname,read.seq))
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
            encompassing_qnames = fields[JUNC_MAP_QNAME_COLUMN].split(SEQ_FIELD_DELIM)
            union_cov = len(spanning_qnames.union(encompassing_qnames))
            intersect_cov = len(spanning_qnames.intersection(encompassing_qnames))
            #both_cov = sum(1 for read_tuple in read_tuples if read_tuple[0] in set(encompassing_qnames))
            fields += [cov, intersect_cov, union_cov, ','.join(map(str,anchors[chimera_id])),
                       SEQ_FIELD_DELIM.join([read_tuple[1] for read_tuple in read_tuples]),
                       SEQ_FIELD_DELIM.join([read_tuple[0] for read_tuple in read_tuples])]
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
