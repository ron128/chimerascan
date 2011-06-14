'''
Created on Jun 3, 2011

@author: mkiyer
'''
import logging
import chimerascan.pysam as pysam


def extend_and_pad_sam(input_fastq_files, input_sam_file, output_sam_file):
    infh = pysam.Samfile(input_sam_file, "r")
    outfh = pysam.Samfile(output_sam_file, "w", template=infh)    

    tagdict = dict(r.tags)
    # TODO: bug in pysam handling CP tag, fix by forcing to integer
    if "CP" in tagdict:
        tagdict["CP"] = int(tagdict["CP"])
    # add additional tags
    tagdict.update(tags)
    r.tags = tagdict.items()

def write_reads_to_bam(reads, bamfh, multihits, keep_unmapped):
    num_unmapped = 0
    num_multihits = 0
    for r in reads:
        if r.is_unmapped:
            xm_tag = r.opt('XM')
            if xm_tag < multihits:
                num_unmapped += 1
                if not keep_unmapped:
                    continue
            else:
                num_multihits += 1
        bamfh.write(r)
    return num_unmapped, num_multihits

def sam_to_bam(input_fastq_file, input_sam_file, output_bam_file, 
               multihits, mode, keep_unmapped=True):
    samfh = pysam.Samfile(input_sam_file, "r")
    if mode == "pe":
        fix_iter = fix_pe_alignment_ordering(samfh, open(input_fastq_file),
                                             is_paired=True)
    elif mode == "pesr":
        fix_iter = fix_pe_sr_alignment_ordering(samfh, open(input_fastq_file))
    elif mode == "sr": 
        fix_iter = fix_pe_alignment_ordering(samfh, open(input_fastq_file),
                                             is_paired=False)
    num_unmapped = 0
    num_multihits = 0
    num_frags = 0
    bamfh = pysam.Samfile(output_bam_file, "wb", template=samfh)
    for frags in fix_iter:
        num_frags += 1
        for reads in frags:
            un, mh = write_reads_to_bam(reads, bamfh, multihits, 
                                        keep_unmapped)
            num_unmapped += un
            num_multihits += mh
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
    parser = OptionParser("usage: %prog [options] --mode <mode> <in1.fq> <in2.fq> <in.sam> <out.bam>")
    parser.add_option("--multihits", type="int", dest="multihits", default=100)
    parser.add_option("--mode", dest="mode", choices=["pe", "pesr", "sr"], default=None)
    parser.add_option("--un", action="store_true", dest="keep_unmapped", default=False)
    options, args = parser.parse_args()
    if options.mode is None:
        parser.error("must specify '--mode' (pe, pesr, sr)")
    input_fastq_file = args[0]
    input_sam_file = args[1]    
    output_bam_file = args[2]
    sam_stdin_to_bam(input_fastq_file,
                     input_sam_file,
                     output_bam_file, 
                     multihits=options.multihits,
                     mode=options.mode,
                     keep_unmapped=options.keep_unmapped)