'''
Created on Jan 22, 2011

@author: mkiyer
'''
import logging
import subprocess
import sys

# local imports
from chimerascan.lib import pysam
from chimerascan.lib.base import get_read_length
from fix_alignment_ordering import fix_alignment_ordering


def align_pe_full(fastq_files, 
                  bowtie_index,
                  output_bam_file, 
                  unaligned_fastq_param,
                  maxmultimap_fastq_param,
                  min_fragment_length=0,
                  max_fragment_length=1000,
                  trim5=0,
                  trim3=0,
                  library_type="fr",
                  num_processors=1, 
                  fastq_format="phred33-quals", 
                  multihits=100, 
                  mismatches=2, 
                  bowtie_bin="bowtie", 
                  bowtie_mode="-n"):
    read_length = get_read_length(fastq_files[0])     
    args = [bowtie_bin, "-q", "-S", 
            "-p", str(num_processors),
            "--%s" % fastq_format,
            "-k", str(multihits),
            "-m", str(multihits),
            bowtie_mode, str(mismatches),
            "--minins", min_fragment_length,
            "--maxins", max_fragment_length,
            "--trim5", trim5,
            "--trim3", trim3,
            "--%s" % library_type,
            "--un", unaligned_fastq_param,
            "--max", maxmultimap_fastq_param]
    # use the entire read length as the "seed" here
    if bowtie_mode == "-n":
        args.extend(["-l", str(read_length)])
    args += [bowtie_index, 
             "-1", fastq_files[0],
             "-2", fastq_files[1]]
    #aligned_sam_file]
    args = map(str, args)
    logging.debug("Bowtie alignment args: %s" % (' '.join(args)))
    aln_p = subprocess.Popen(args, stdout=subprocess.PIPE)
    # pipe the bowtie SAM output to a filter that writes BAM format
    args = [sys.executable, __file__, output_bam_file, fastq_files[0], str(multihits)]
    logging.debug("SAM to BAM converter args: %s" % (' '.join(args)))
    retcode = subprocess.call(args, stdin=aln_p.stdout)
    if retcode != 0:
        return retcode
    return aln_p.wait()

def sam_stdin_to_bam(output_bam_file, input_fastq_file, multihits):
    samfh = pysam.Samfile("-", "r")
    bamfh = pysam.Samfile(output_bam_file, "wb", template=samfh)
    
    for pe_reads in fix_alignment_ordering(samfh, open(input_fastq_file), 
                                           is_paired=True):
        for reads in pe_reads:
            for r in reads:
                bamfh.write(r)
#        if r.is_unmapped:
#            xm_tag = r.opt('XM')
#            # keep multihits in the BAM file but remove nonmapping reads
#            # since these will specifically be remapped later
#            if xm_tag < multihits:
#                num_unmapped += 1
#                continue
#            num_multihits += 1
    bamfh.close()
    samfh.close()
    #logging.debug("[SAMTOBAM] Filtered %d unmapped reads" % (num_unmapped))
    #logging.debug("[SAMTOBAM] Allowed %d highly multimapping reads to pass through as unmapped" % (num_multihits))
    logging.info("[SAMTOBAM] Finished converting SAM -> BAM")

if __name__ == '__main__':
    sam_stdin_to_bam(sys.argv[1], sys.argv[2], int(sys.argv[3]))
