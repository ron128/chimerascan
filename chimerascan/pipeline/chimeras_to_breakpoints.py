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

DEFAULT_HOMOLOGY_MISMATCHES = config.BREAKPOINT_HOMOLOGY_MISMATCHES

def calc_homology(seq1, seq2, num_mismatches):
    smallest_len = min(len(seq1), len(seq2))
    mm = 0
    i = 0
    for i in xrange(smallest_len):
        if seq1[i] != seq2[i]:
            mm += 1
            if mm > num_mismatches:
                break
    return i

def determine_chimera_breakpoints(index_dir, read_length, 
                                  input_chimera_file, output_chimera_file, 
                                  breakpoint_map_file, breakpoint_fasta_file,
                                  homology_mismatches=DEFAULT_HOMOLOGY_MISMATCHES):
    # open the reference sequence fasta file
    ref_fasta_file = os.path.join(index_dir, config.ALIGN_INDEX + ".fa")
    ref_fa = pysam.Fastafile(ref_fasta_file)
    # output files
    chimerafh = open(output_chimera_file, "w")
    breakpointfh = open(breakpoint_map_file, "w")
    fasta_output_fh = open(breakpoint_fasta_file, "w")
    breakpoints = collections.defaultdict(lambda: [])
    breaknum = 0
    for c in Chimera.parse(open(input_chimera_file)):
        # retrieve transcript coordinates of 5' and 3' partners
        ref5p = config.GENE_REF_PREFIX + c.partner5p.tx_name
        ref3p = config.GENE_REF_PREFIX + c.partner3p.tx_name
        start5p, end5p = c.partner5p.start, c.partner5p.end
        start3p, end3p = c.partner3p.start, c.partner3p.end
        # get intervals for breakpoint sequence
        breakpoint_start5p = max(start5p, end5p - read_length + 1)
        breakpoint_end3p = min(end3p, start3p + read_length - 1)
        # fetch sequence
        seq5p = ref_fa.fetch(ref5p, breakpoint_start5p, end5p)
        seq3p = ref_fa.fetch(ref3p, start3p, breakpoint_end3p)
        if len(seq5p) < read_length - 1:
            logging.warning("Could not extract sequence of length >%d from "
                            "5' partner of chimera %s, only retrieved "
                            "sequence of %d" % 
                            (read_length-1, c.name, len(seq5p)))
            # pad sequence
            padding = (read_length - 1) - len(seq5p)
            seq5p = ("N" * padding) + seq5p
        if len(seq3p) < read_length - 1:
            logging.warning("Could not extract sequence of length >%d from "
                            "3' partner of chimera %s, only retrieved "
                            "sequence of %d" % 
                            (read_length-1, c.name, len(seq3p)))
            # pad sequence
            padding = (read_length - 1) - len(seq3p)
            seq3p = seq3p + ("N" * padding)
        # fetch continuation sequence of non-fusion gene
        homolog_end5p = end5p + read_length - 1
        homolog_start3p = max(0, start3p - read_length + 1)
        homolog5p = ref_fa.fetch(ref3p, homolog_start3p, start3p)
        homolog3p = ref_fa.fetch(ref5p, end5p, homolog_end5p)
        # find homology between 5' gene and 3' gene
        homology_length_5p = calc_homology(seq5p[::-1], homolog5p[::-1], 
                                           homology_mismatches)
        homology_length_3p = calc_homology(seq3p, homolog3p, 
                                           homology_mismatches)        
        # create a Breakpoint and add to dictionary
        seq = seq5p + seq3p
        if seq in breakpoints:
            b = breakpoints[seq]
        else:
            b = Breakpoint()
            b.name = "B%07d" % (breaknum)
            breaknum += 1
            b.seq5p = seq5p
            b.seq3p = seq3p
            breakpoints[seq] = b
        # add sequence to dictionary and group fusion candidates together
        # if they have the same location and junction sequence
        b.chimera_names.append(c.name)
        # update Chimera object with breakpoint information
        c.breakpoint_name = b.name
        c.breakpoint_homology_5p = homology_length_5p
        c.breakpoint_homology_3p = homology_length_3p
        # write Chimera
        fields = c.to_list()
        print >>chimerafh, '\t'.join(map(str, c.to_list()))
    # now extract the unique junction sequences
    # and write them to a fasta file
    for seq,b in breakpoints.iteritems():
        # write to fasta file
        print >>fasta_output_fh, ">%s\n%s" % (b.name, seq)
        # write to breakpoint map file
        fields = b.to_list()
        print >>breakpointfh, '\t'.join(map(str, fields))
    # close files
    fasta_output_fh.close()
    breakpointfh.close()
    chimerafh.close()

def main():
    from optparse import OptionParser
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = OptionParser("usage: %prog [options] <index> <read_length> "
                          "<chimeras.txt> <chimeras.out.txt> " 
                          "<breakpoints.txt> <breakpoints.fa>")
    parser.add_option("--homology-mismatches", type="int", 
                      dest="homology_mismatches", 
                      default=config.BREAKPOINT_HOMOLOGY_MISMATCHES,
                      help="Number of mismatches to tolerate when computing "
                      "homology between gene and its chimeric partner "
                      "[default=%default]")
    options, args = parser.parse_args()
    index_dir = args[0]
    read_length = int(args[1])
    input_chimera_file = args[2]
    output_chimera_file = args[3]
    breakpoint_map_file = args[4]
    breakpoint_fasta_file = args[5]
    determine_chimera_breakpoints(index_dir, 
                                  read_length, 
                                  input_chimera_file, 
                                  output_chimera_file, 
                                  breakpoint_map_file, 
                                  breakpoint_fasta_file,
                                  homology_mismatches=options.homology_mismatches)

if __name__ == '__main__':
    main()
