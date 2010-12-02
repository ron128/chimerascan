'''
Created on Nov 7, 2010

@author: mkiyer
'''
import logging
import argparse
import collections
import numpy as np
from scipy.stats import chisquare

def read_chimera_mapping_file(filename):
    chimera_refs = collections.defaultdict(lambda: [])
    for line in open(filename):
        fields = line.strip().split('\t')
        chimera_id = fields[0]
        bedpe_fields = fields[1:]        
        chimera_refs[chimera_id].append(bedpe_fields)
    return dict(chimera_refs)

class BowtieAlignedRead(object):
    __slots__ = ('qname', 'mate', 'strand', 'rname', 'pos', 'seq', 'qual', 'md')
    @staticmethod
    def from_bowtie_line(line):
        fields = line.strip().split('\t')
        read = BowtieAlignedRead()
        read.qname = fields[0][:-2]
        read.mate = int(fields[0][-1])
        read.strand = fields[1]
        read.rname = fields[2]
        read.pos = int(fields[3])
        read.seq = fields[4]
        read.qual = fields[5]
        read.md = []
        if len(fields) == 8:
            mismatch_desc = fields[7]
            for mismatch_string in mismatch_desc.split(','):
                pos,basechange = mismatch_string.split(':')
                from_base,to_base = basechange.split('>')
                read.md.append((int(pos), from_base, to_base))
        return read

def parse_bowtie_multihits(bowtie_output_file):
    reads = []
    for line in open(bowtie_output_file):
        read = BowtieAlignedRead.from_bowtie_line(line)
        if len(reads) != 0 and reads[-1].qname != read.qname:
            yield reads
            reads = []
        reads.append(read)
    if len(reads) != 0:
        yield reads
    
def join_spanning_reads(bowtie_output_files, chimera_refs,
                        original_read_length,
                        anchor_min, anchor_max,
                        max_anchor_mismatches):
    chimera_counts = collections.defaultdict(lambda: 0)
    chimera_reads = collections.defaultdict(lambda: [])
    anchors = {}
    junc_positions = {}
    max_anchor = int((original_read_length + 1)/2)
    for chimera_id in chimera_refs.iterkeys():
        left_junc_length = int(chimera_id.split(":")[1])
        junc_positions[chimera_id] = left_junc_length
        anchors[chimera_id] = np.zeros(max_anchor, dtype=np.int)
    # count read alignments to chimera junctions
    for bowtie_output_file in bowtie_output_files:
        for reads in parse_bowtie_multihits(bowtie_output_file):
            filtered_reads = []
            for read in reads:
                # filter out reads with small anchor sequence
                read_length = len(read.seq)                
                junc_pos = junc_positions[read.rname]
                left_anchor_bp = junc_pos - read.pos
                right_anchor_bp = read.pos + read_length - junc_pos
                # keep track of minimum number of bp spanning the chimera
                anchor = min(left_anchor_bp, right_anchor_bp)
                anchors[read.rname][anchor - 1] += 1
                # see if read passes filters 
                read_ok = True
                if anchor < anchor_min:
                    read_ok = False
                elif anchor < anchor_max:
                    # check mismatch descriptors
                    if left_anchor_bp < anchor_max:
                        anchor_interval = (0, left_anchor_bp)
                    else:
                        anchor_interval = (read_length - right_anchor_bp, read_length)                        
                    mismatches = [mm for mm in read.md
                                  if mm[0] >= anchor_interval[0] and mm[0] < anchor_interval[1]]
                    if len(mismatches) > max_anchor_mismatches:
                        read_ok = False
                if read_ok:
                    filtered_reads.append(read)
                #else:
                #    print read.qname, read.mate, read.strand, read.rname, read.pos, read.md
            # count reads that pass anchor filter
            if len(filtered_reads) == 0:
                continue
            #weight = 1.0 / len(filtered_reads)
            for read in filtered_reads:
                chimera_counts[read.rname] += 1
                chimera_reads[read.rname].append((read.qname,read.seq))
    # assign junction coverage to chimera candidates
    for chimera_id,bedpe_records in chimera_refs.iteritems():
        cov = chimera_counts[chimera_id]
        read_tuples = chimera_reads[chimera_id]
        if len(read_tuples) == 0:
            read_tuples = [("None", "None")]
        for fields in bedpe_records:
            # intersect encompassing with spanning reads to see overlap
            encompassing_qnames = fields[16].split(',')[:-1]
            #both_qnames = set(qnames).intersection(encompassing_qnames)
            both_cov = sum(1 for read_tuple in read_tuples if read_tuple[0] in set(encompassing_qnames))
            #print '\t'.join(map(str, fields))
            #print observed_arr
            #print expected_arr
            #print cov, csq, pval
            #continue
            fields += [cov, both_cov, ','.join(map(str,anchors[chimera_id])),
                       ','.join(['|'.join(read_tuple) for read_tuple in read_tuples])] 
#                       ','.join(qnames)]
            yield fields

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")    
    parser = argparse.ArgumentParser()
    parser.add_argument("--rlen", type=int, dest="read_length")
    parser.add_argument("--anchor-min", type=int, dest="anchor_min", default=0)
    parser.add_argument("--anchor-max", type=int, dest="anchor_max", default=0)
    parser.add_argument("--anchor-mismatches", type=int, dest="anchor_mismatches", default=0)
    parser.add_argument("chimera_mapping_file")
    parser.add_argument("spanning_chimera_file")
    parser.add_argument("bowtie_files", nargs="+")
    options = parser.parse_args()
    
    chimera_refs = read_chimera_mapping_file(options.chimera_mapping_file)
    f = open(options.spanning_chimera_file, "w")
    for bedpe_fields in join_spanning_reads(options.bowtie_files, chimera_refs,
                                            options.read_length,
                                            options.anchor_min,
                                            options.anchor_max, 
                                            options.anchor_mismatches):
        print >>f, '\t'.join(map(str, bedpe_fields))
    f.close()

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
