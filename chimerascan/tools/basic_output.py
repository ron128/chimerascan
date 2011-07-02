'''
Created on Jul 1, 2011

@author: mkiyer
'''
import logging

from chimerascan.lib.chimera import Chimera

def main():
    from optparse import OptionParser
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = OptionParser("usage: %prog [options] <in.txt>")
    options, args = parser.parse_args()
    input_file = args[0]    
    for c in Chimera.parse(open(input_file)):
        fields = [c.partner5p.tx_name, c.partner3p.tx_name,
                  '|'.join([c.partner5p.gene_name, c.partner3p.gene_name]),
                  c.get_weighted_cov(), c.get_unique_spanning_reads()]
        print '\t'.join(map(str, fields))


if __name__ == '__main__':
    main()
