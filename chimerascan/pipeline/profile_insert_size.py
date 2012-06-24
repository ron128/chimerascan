'''
Created on Jan 24, 2011

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
import sys
import pysam

# local imports
from chimerascan.lib.fragment_size_distribution import InsertSizeDistribution

def main():
    from optparse import OptionParser
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = OptionParser("usage: %prog [options] <bam> <out.bedpe>")
    parser.add_option('--min-fragment-length', dest="min_fragment_length", 
                      type="int", default=0)
    parser.add_option('--max-fragment-length', dest="max_fragment_length", 
                      type="int", default=1000)
    parser.add_option('--max-samples', dest="max_samples", 
                      type="int", default=None)
    parser.add_option('-o', dest="output_file", default=None) 
    options, args = parser.parse_args()
    input_bam_file = args[0]
    bamfh = pysam.Samfile(input_bam_file, "rb")
    isizedist = InsertSizeDistribution.from_bam(bamfh, options.min_fragment_length, 
                                                options.max_fragment_length, 
                                                options.max_samples)
    bamfh.close()
    if options.output_file is not None:
        f = open(options.output_file, "w")
    else:
        f = sys.stdout
    isizedist.to_file(f)
    if options.output_file is not None:
        f.close()
    logging.info("Insert size samples=%d mean=%f std=%f median=%d mode=%d" % 
                 (isizedist.n, isizedist.mean(), isizedist.std(), 
                  isizedist.percentile(50.0), isizedist.mode()))
    

if __name__ == '__main__':
    main()