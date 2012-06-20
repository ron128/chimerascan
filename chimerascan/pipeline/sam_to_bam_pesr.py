'''
Created on Jun 1, 2012

@author: mkiyer
'''
'''
Created on Jun 2, 2011

@author: mkiyer
'''
import sys
import logging
import argparse

# local imports
import chimerascan.pysam as pysam
from chimerascan.lib import config
from chimerascan.lib.sam import soft_pad_read
from chimerascan.lib.seq import parse_fastq_record

def sam_to_bam_pesr(input_sam_file, input_fastq_file, output_bam_file):
    samfh = pysam.Samfile(input_sam_file, "r")
    fqiter = parse_fastq_record(open(input_fastq_file))
    bamfh = pysam.Samfile(output_bam_file, "wb", template=samfh)
    num_reads = 0
    # get first fastq record
    fqrec = fqiter.next()
    for r in samfh:
        # get next fastq record
        if r.qname != fqrec.qname:
            fqrec = fqiter.next()
        # TODO: eventually remove assert statement
        assert r.qname == fqrec.qname
        # remove mate number from qname
        r.qname = r.qname[1:]
        # pad read that has been trimmed
        soft_pad_read(fqrec, r)
        # reset paired-end flags for read
        r.is_paired = True
        if fqrec.readnum == 1:
            r.is_read1 = True
        elif fqrec.readnum == 2:
            r.is_read2 = True
        else:
            assert False
        # TODO: remove is_secondary bit here because this does not 
        # make sense in the context of paired-end alignments
        r.is_secondary = False
        bamfh.write(r)
        num_reads += 1
    logging.debug("Found %d reads" % (num_reads))
    bamfh.close()
    samfh.close()
    return config.JOB_SUCCESS

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("input_sam_file")
    parser.add_argument("input_fastq_file")
    parser.add_argument("output_bam_file") 
    args = parser.parse_args()
    return sam_to_bam_pesr(args.input_sam_file, args.input_fastq_file, args.output_bam_file)

if __name__ == '__main__':
    sys.exit(main())
