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
    num_frags = 0
    # get first fastq record
    fqrec = fqiter.next()
    for r in samfh:
        # get next fastq record
        if r.qname != fqrec.qname:
            fqrec = fqiter.next()
        # TODO: eventually remove assert statement
        assert r.qname == fqrec.qname
        # reset paired-end flags for read
        r.is_paired = True
        if fqrec.readnum == 1:
            r.is_read1 = True
        elif fqrec.readnum == 2:
            r.is_read2 = True
        else:
            assert False
        r.qname = r.qname[1:]
        # pad read that has been trimmed
        soft_pad_read(fqrec, r)
        bamfh.write(r)
        num_frags += 1
    logging.debug("Found %d fragments" % (num_frags))
    bamfh.close()
    samfh.close()
    return config.JOB_SUCCESS

#def sam_to_bam(input_fastq_files, input_sam_file, output_bam_file, 
#               quals, multihits, pe_sr_mode=False, softclip=True, 
#               keep_unmapped=True):
#    samfh = pysam.Samfile(input_sam_file, "r")
#    num_unmapped = 0
#    num_multihits = 0
#    num_frags = 0
#    bamfh = pysam.Samfile(output_bam_file, "wb", template=samfh)
#    # setup fastq parsing
#    if softclip and (quals != SANGER_FORMAT):
#        kwargs = {"convert_quals": True, "qual_format": quals}
#    else:
#        kwargs = {"convert_quals": False}
#    fqiters = [parse_fastq_record(open(fq), **kwargs) for fq in input_fastq_files]
#    
#    # handle single-read and paired-end
#    if len(fqiters) == 1:
#        reorder_func = fix_sr_alignment_ordering(samfh, fqiters[0])
#    else:
#        reorder_func = fix_alignment_ordering(samfh, fqiters, pe_sr_mode)
#    # iterate through buffer
#    for bufitems in reorder_func:
#        num_frags += 1
#        for bufitem in bufitems:
#            for r in bufitem.reads:
#                # softclip uses the fastq record to replace the sequence
#                # and quality scores of the read 
#                if softclip:
#                    soft_pad_read(bufitem.fqrec, r)
#                # keep statistics of unmapped/multimapped reads and
#                # suppress output if 'keep_unmapped' is False
#                if r.is_unmapped:
#                    xm_tag = r.opt('XM')
#                    if xm_tag < multihits:
#                        num_unmapped += 1
#                        if not keep_unmapped:
#                            continue
#                    else:
#                        num_multihits += 1
#                bamfh.write(r)
#    for fqfh in fqiters:
#        fqfh.close()
#    bamfh.close()
#    samfh.close()
#    logging.debug("Found %d fragments" % (num_frags))
#    logging.debug("\t%d unmapped reads" % (num_unmapped))
#    logging.debug("\t%d multimapping (>%dX) reads" % 
#                  (num_multihits, multihits))

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
