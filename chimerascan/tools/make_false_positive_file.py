'''
Created on Jul 6, 2011

@author: mkiyer
'''
import logging
import sys
import collections

from chimerascan.lib.chimera import Chimera

def main():
    from optparse import OptionParser
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = OptionParser("usage: %prog [options] <chimeras.txt> [<chimeras2.txt> <chimeras3.txt> ...]")
    parser.add_option("-o", dest="output_file", default=None,
                      help="output file [default=stdout]")
    parser.add_option("-n", dest="num_files", type="int", default=1,
                      help="chimera must be recurrent in N samples "
                      "to make considered a false positive "
                      "[default=%default]")
    options, args = parser.parse_args()
    input_files = args
    false_pos_chimeras = collections.defaultdict(lambda: 0)
    for input_file in input_files:
        logging.info("Processing file %s" % (input_file))
        num_chimeras = 0
        for c in Chimera.parse(open(input_file)):
            key = (c.partner5p.tx_name, c.partner5p.end, c.partner3p.tx_name, c.partner3p.start)
            false_pos_chimeras[key] += 1
            num_chimeras += 1
        logging.info("\tchimeras in file: %d" % (num_chimeras))
        logging.info("\tcurrent false positive candidates: %d" % (len(false_pos_chimeras)))
    if options.output_file is None:
        fileh = sys.stdout
    else:
        fileh = open(options.output_file, "w")
    for key,recurrence in false_pos_chimeras.iteritems():
        if recurrence >= options.num_files:
            print >>fileh, '\t'.join(map(str,key))
    if options.output_file is not None:
        fileh.close()

if __name__ == '__main__':
    main()