'''
Created on Jul 14, 2011

@author: mkiyer
'''
import logging
import sys
import random

def split_seq(seq, chars_per_line):
    pos = 0
    newseq = []
    while pos < len(seq):
        if pos + chars_per_line > len(seq):        
            endpos = len(seq)
        else:
            endpos = pos + chars_per_line
        newseq.append(seq[pos:endpos])
        pos = endpos
    return '\n'.join(newseq)

def main():
    from optparse import OptionParser
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = OptionParser("usage: %prog [options] <fasta_id>")
    parser.add_option("-o", dest="output_file", default=None,
                      help="output file [default=stdout]")
    parser.add_option("-l", dest="length", type="int", default=int(1e6),
                      help="number of bases of random dna to generate "
                      "[default=%default]")
    options, args = parser.parse_args()
    if options.output_file is None:
        fileh = sys.stdout
    else:
        fileh = open(options.output_file, "w")
    # make a string of dna (not computationally efficient)
    dna_lookup = ["A", "T", "G", "C"]
    dna = ''.join([dna_lookup[random.randrange(0,4)] 
                   for x in xrange(options.length)])   
    print >>fileh, ">%s\n%s" % (args[0], split_seq(dna,50))
    if options.output_file is not None:
        fileh.close()


if __name__ == '__main__':
    main()