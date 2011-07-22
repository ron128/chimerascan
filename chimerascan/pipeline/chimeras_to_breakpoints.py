'''
Created on Jun 11, 2011

@author: mkiyer
'''
import logging
import os
import collections

from chimerascan import pysam
from chimerascan.lib import config
from chimerascan.lib.breakpoint import Breakpoint
from chimerascan.lib.chimera import Chimera
from chimerascan.lib.batch_sort import batch_sort
from chimerascan.lib.seq import split_seq

def chimeras_to_breakpoints(input_file, breakpoint_map_file, breakpoint_fasta_file,
                            tmp_dir):
    # TODO: use make_temp
    # write breakpoint/chimera relationships to file
    tmp_file = os.path.join(tmp_dir, "breakpoint_info.tmp")
    f = open(tmp_file, "w")
    for c in Chimera.parse(open(input_file)):
        # write to temp file
        print >>f, "%s\t%s" % (c.breakpoint_seq_5p + c.breakpoint_seq_3p, c.name)
    f.close()
    # sort breakpoint file
    def sortfunc(line):
        fields = line.strip().split('\t')
        return fields[0]
    tempdirs = [tmp_dir]
    sorted_tmp_file = os.path.join(tmp_dir, "breakpoint_info.srt.tmp")
    batch_sort(input=tmp_file,
               output=sorted_tmp_file,
               key=sortfunc,
               buffer_size=32000,
               tempdirs=tempdirs)
    os.remove(tmp_file)
    # parse and build breakpoint -> chimera map
    fastafh = open(breakpoint_fasta_file, "w")
    mapfh = open(breakpoint_map_file, "w")
    breakpoint_num = 0
    prev_seq = None
    chimera_names = set()
    for line in open(sorted_tmp_file):
        fields = line.strip().split('\t')
        seq = fields[0]
        chimera_name = fields[1]
        if seq != prev_seq:
            if len(chimera_names) > 0:
                # write to fasta
                name = "B%07d" % (breakpoint_num)
                print >>fastafh, ">%s\n%s" % (name, split_seq(prev_seq))
                print >>mapfh, "%s\t%s\t%s" % (name, prev_seq, ",".join(chimera_names))
                chimera_names = set()
                breakpoint_num += 1
            prev_seq = seq
        chimera_names.add(chimera_name) 
    if len(chimera_names) > 0:
        name = "B%07d" % (breakpoint_num)
        print >>fastafh, ">%s\n%s" % (name, split_seq(prev_seq))
        print >>mapfh, "%s\t%s\t%s" % (name, prev_seq, ",".join(chimera_names))
    os.remove(sorted_tmp_file)


def main():
    from optparse import OptionParser
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = OptionParser("usage: %prog [options] <chimeras.bedpe> "
                          "<breakpoints.txt> <breakpoints.fa> <tmp_dir>") 
    options, args = parser.parse_args()
    input_file = args[0]
    breakpoint_map_file = args[1]
    breakpoint_fasta_file = args[2]
    tmp_dir = args[3]
    chimeras_to_breakpoints(input_file, breakpoint_map_file, breakpoint_fasta_file, tmp_dir)


if __name__ == '__main__':
    main()
