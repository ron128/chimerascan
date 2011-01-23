'''
Created on Jan 23, 2011

@author: mkiyer
'''
import logging
import re

# local imports
import pysam
from find_discordant_reads import Chimera

def parse_qname_file(qname_fh):
    for line in qname_fh:
        qname, mate = line.strip().split('\t')
        yield qname, int(mate)

def parse_discordant_reads(infh):
    for line in infh:
        chimera = Chimera.from_bedpe(line)
        yield chimera

def parse_fastq(line_iter):
    try:        
        qname = line_iter.next().rstrip()[1:]
        newqname = re.split(r'/\d$', qname)[0]
        suffix_length = len(qname) - len(newqname)                    
        seq = line_iter.next().rstrip()
        line_iter.next()
        qual = line_iter.next().rstrip()
        yield newqname, seq, qual
        while True:
            # qname
            qname = line_iter.next().rstrip()[1:]
            qname = qname[:len(qname)-suffix_length]
            # seq
            seq = line_iter.next().rstrip()
            # qname again (skip)
            line_iter.next()
            # qual
            qual = line_iter.next().rstrip()
            yield qname, seq, qual
    except StopIteration:
        pass

def sam_stdin_to_bam(output_bam_file, multihits):
    samfh = pysam.Samfile("-", "r")
    bamfh = pysam.Samfile(output_bam_file, "wb", template=samfh)
    num_unmapped = 0
    num_multihits = 0
    for r in samfh:
        if r.is_unmapped:
            xm_tag = r.opt('XM')
            # keep multihits in the BAM file but remove nonmapping reads
            # since these will specifically be remapped later
            if xm_tag < multihits:
                num_unmapped += 1
                continue
            num_multihits += 1
        bamfh.write(r)
    bamfh.close()
    samfh.close()
    logging.debug("[SAMTOBAM] Filtered %d unmapped reads" % (num_unmapped))
    logging.debug("[SAMTOBAM] Allowed %d highly multimapping reads to pass through as unmapped" % (num_multihits))
    logging.info("[SAMTOBAM] Finished converting SAM -> BAM")



def extend_sequences(input_fastq_files, discordant_bedpe_file,
                     output_discordant_bedpe_file, 
                     output_spanning_fastq_file):
    fastq_iters = [parse_fastq(open(fq)) for fq in input_fastq_files]
    discordant_iter = parse_discordant_reads(open(discordant_bedpe_file))     
    bedpe_fh = open(output_discordant_bedpe_file, "w")
    fastq_fh = open(output_spanning_fastq_file, "w")    
    try:
        
        pass
    except StopIteration:
        pass
    
    

def main():
    from optparse import OptionParser
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = OptionParser("usage: %prog [options]")
    options, args = parser.parse_args()    
    discordant_bedpe_file = args[0]
    output_discordant_bedpe_file = args[1]
    output_spanning_fastq_file = args[2]
    input_fastq_files = args[3:]
    extend_sequences(input_fastq_files, discordant_bedpe_file,
                     output_discordant_bedpe_file, 
                     output_spanning_fastq_file)

if __name__ == '__main__':
    main()

