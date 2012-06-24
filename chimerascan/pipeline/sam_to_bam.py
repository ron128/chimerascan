'''
Created on Jun 2, 2011

@author: mkiyer
'''
import sys
import logging
import argparse

import pysam

# local imports
from chimerascan.lib import config

def sam_to_bam(input_sam_file, output_bam_file):
    samfh = pysam.Samfile(input_sam_file, "r")
    bamfh = pysam.Samfile(output_bam_file, "wb", template=samfh)
    num_frags = 0
    for r in samfh:
        bamfh.write(r)
        num_frags += 1
    logging.debug("Found %d fragments" % (num_frags))
    samfh.close()
    bamfh.close()
    return config.JOB_SUCCESS

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("input_sam_file")
    parser.add_argument("output_bam_file") 
    args = parser.parse_args()
    return sam_to_bam(args.input_sam_file, args.output_bam_file)

if __name__ == '__main__':
    sys.exit(main())
