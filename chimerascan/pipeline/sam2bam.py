'''
Created on Jun 2, 2011

@author: mkiyer
'''
import logging

# local imports
import chimerascan.pysam as pysam
from chimerascan.lib.fix_alignment_ordering import fix_alignment_ordering, fix_sr_alignment_ordering
from chimerascan.lib.sam import soft_pad_read
from chimerascan.lib.seq import FASTQ_QUAL_FORMATS, SANGER_FORMAT, parse_fastq_record

def sam_to_bam(input_fastq_files, input_sam_file, output_bam_file, 
               quals, multihits, pe_sr_mode=False, softclip=True, 
               keep_unmapped=True):
    samfh = pysam.Samfile(input_sam_file, "r")
    num_unmapped = 0
    num_multihits = 0
    num_frags = 0
    bamfh = pysam.Samfile(output_bam_file, "wb", template=samfh)
    # setup fastq parsing
    if softclip and (quals != SANGER_FORMAT):
        kwargs = {"convert_quals": True, "qual_format": quals}
    else:
        kwargs = {"convert_quals": False}
    fqiters = [parse_fastq_record(open(fq), **kwargs) for fq in input_fastq_files]
    
    # handle single-read and paired-end
    if len(fqiters) == 1:
        reorder_func = fix_sr_alignment_ordering(samfh, fqiters[0])
    else:
        reorder_func = fix_alignment_ordering(samfh, fqiters, pe_sr_mode)
    # iterate through buffer
    for bufitems in reorder_func:
        num_frags += 1
        for bufitem in bufitems:
            for r in bufitem.reads:
                # softclip uses the fastq record to replace the sequence
                # and quality scores of the read 
                if softclip:
                    soft_pad_read(bufitem.fqrec, r)
                # keep statistics of unmapped/multimapped reads and
                # suppress output if 'keep_unmapped' is False
                if r.is_unmapped:
                    xm_tag = r.opt('XM')
                    if xm_tag < multihits:
                        num_unmapped += 1
                        if not keep_unmapped:
                            continue
                    else:
                        num_multihits += 1
                bamfh.write(r)
    for fqfh in fqiters:
        fqfh.close()
    bamfh.close()
    samfh.close()
    logging.debug("Found %d fragments" % (num_frags))
    logging.debug("\t%d unmapped reads" % (num_unmapped))
    logging.debug("\t%d multimapping (>%dX) reads" % 
                  (num_multihits, multihits))

if __name__ == '__main__':
    from optparse import OptionParser
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = OptionParser("usage: %prog [options] <out.bam> <in.sam> <in1.fq> [<in2.fq>]")
    parser.add_option("--multihits", type="int", dest="multihits", default=100)
    parser.add_option("--quals", dest="quals", 
                      choices=FASTQ_QUAL_FORMATS,
                      default=SANGER_FORMAT)
    parser.add_option("--pesr", action="store_true", dest="pe_sr_mode", default=False)
    parser.add_option("--softclip", action="store_true", dest="softclip", default=False)
    parser.add_option("--un", action="store_true", dest="keep_unmapped", default=False)
    options, args = parser.parse_args()
    output_bam_file = args[0]
    input_sam_file = args[1]    
    input_fastq_files = args[2:]
    sam_to_bam(input_fastq_files,
               input_sam_file,
               output_bam_file, 
               quals=options.quals,
               multihits=options.multihits,
               pe_sr_mode=options.pe_sr_mode,
               softclip=options.softclip,
               keep_unmapped=options.keep_unmapped)
