'''
Created on Jun 11, 2011

@author: mkiyer
'''
import logging

from chimerascan.lib.chimera import Chimera
from chimerascan.lib.batch_sort import batch_sort
from chimerascan.lib.seq import split_seq

def chimeras_to_breakpoints(input_file, breakpoint_sorted_chimera_file, 
                            breakpoint_map_file, breakpoint_fasta_file,
                            tmp_dir):
    # sort chimera file by breakpoint name
    def sortfunc(line):
        fields = line.strip().split('\t')
        return fields[Chimera.BREAKPOINT_NAME_FIELD]
    tempdirs = [tmp_dir]
    batch_sort(input=input_file,
               output=breakpoint_sorted_chimera_file,
               key=sortfunc,
               buffer_size=32000,
               tempdirs=tempdirs)
    # parse and build breakpoint -> chimera map
    fastafh = open(breakpoint_fasta_file, "w")
    mapfh = open(breakpoint_map_file, "w")
    prev_breakpoint_name = None
    prev_seq = None
    chimera_names = set()
    for c in Chimera.parse(open(breakpoint_sorted_chimera_file)):        
        seq = c.breakpoint_seq_5p + c.breakpoint_seq_3p
        if c.breakpoint_name != prev_breakpoint_name:
            if len(chimera_names) > 0:
                # write to fasta
                print >>fastafh, ">%s\n%s" % (prev_breakpoint_name, split_seq(prev_seq))
                # write to map file
                print >>mapfh, "%s\t%s\t%s" % (prev_breakpoint_name, 
                                               prev_seq, 
                                               ",".join(sorted(chimera_names)))
                chimera_names = set()
            prev_seq = seq
            prev_breakpoint_name = c.breakpoint_name
        chimera_names.add(c.name)
    if len(chimera_names) > 0:
        print >>fastafh, ">%s\n%s" % (prev_breakpoint_name, split_seq(prev_seq))
        print >>mapfh, "%s\t%s\t%s" % (prev_breakpoint_name, prev_seq, ",".join(chimera_names))
    fastafh.close()
    mapfh.close()


def main():
    from optparse import OptionParser
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = OptionParser("usage: %prog [options] <chimeras.bedpe> "
                          "<sorted_chimeras.bedpe> "
                          "<breakpoints.txt> <breakpoints.fa> <tmp_dir>") 
    options, args = parser.parse_args()
    input_file = args[0]
    breakpoint_sorted_chimera_file = args[1]
    breakpoint_map_file = args[2]
    breakpoint_fasta_file = args[3]
    tmp_dir = args[3]
    chimeras_to_breakpoints(input_file, breakpoint_sorted_chimera_file,
                            breakpoint_map_file, breakpoint_fasta_file, tmp_dir)


if __name__ == '__main__':
    main()
