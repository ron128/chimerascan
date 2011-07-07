'''
Created on Jul 6, 2011

@author: mkiyer
'''
import logging
import sys

from chimerascan.lib.chimera import Chimera

def main():
    from optparse import OptionParser
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = OptionParser("usage: %prog [options] <chimeras.txt> [<chimeras2.txt> <chimeras3.txt> ...]")
    parser.add_option("-o", dest="output_file", default=None,
                      help="output file [default=stdout]")
    options, args = parser.parse_args()
    input_files = args
    false_pos_pairs = set()
    if options.output_file is None:
        fileh = sys.stdout
    else:
        fileh = open(options.output_file, "w")
    for input_file in input_files:
        logging.info("Processing file %s" % (input_file))
        num_chimeras = 0
        for c in Chimera.parse(open(input_file)):
            key = (c.partner5p.tx_name, c.partner5p.end, c.partner3p.tx_name, c.partner3p.start)
            if key not in false_pos_pairs:                
                print >>fileh, '\t'.join(map(str,key))
                false_pos_pairs.add(key)
            num_chimeras += 1
        logging.info("\tchimeras in file: %d" % (num_chimeras))
        logging.info("\tcurrent false positive candidates: %d" % (len(false_pos_pairs)))
    if options.output_file is not None:
        fileh.close()

if __name__ == '__main__':
    main()