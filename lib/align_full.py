'''
Created on Jan 22, 2011

@author: mkiyer
'''
import logging
import subprocess
import sys

# local imports
import pysam
from base import get_read_length

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
    args = [sys.executable, __file__, output_bam_file, str(multihits)]
    logging.debug("SAM to BAM converter args: %s" % (' '.join(args)))
    retcode = subprocess.call(args, stdin=aln_p.stdout)
    if retcode != 0:
        return retcode
    return aln_p.wait()

def fix_alignment_order(samfh, is_paired=True, maxlen=100000):
    num_mates = 2 if is_paired else 1    
    # function for initializing new buffer list
    buf_init_func = lambda: tuple(list() for m in xrange(num_mates))
    qname_ind_map = {}
    start_ind = 0
    end_ind = 0
    cur_ind = 0
    buf = [None] * maxlen
    for read in samfh:
        qname = read.qname
        mate = 0 if read.is_read1 else 1
        # see whether this read is already in the buffer
        if qname not in qname_ind_map:
            buf_size = len(qname_ind_map)
            # test if buffer has become full
            if buf_size == maxlen:
                # buffer full so return first read
                yield buf[start_ind]
                # delete the read qname to decrease the buffer size by
                # one and allow parser to iterate until it is full again
                return_qname = buf[start_ind][0][0].qname
                del qname_ind_map[return_qname]
                # advance start index
                start_ind += 1
                if start_ind == maxlen:
                    start_ind = 0
            # reset end index in buffer            
            qname_ind_map[qname] = end_ind
            buf[end_ind] = buf_init_func()
            # get current index for insertion
            cur_ind = end_ind
            # advance end index
            end_ind += 1
            if end_ind == maxlen:
                end_ind = 0
        else:
            # grab buffer index for this read qname
            cur_ind = qname_ind_map[qname]
        # add read to buffer
        buf[cur_ind][mate].append(read)        
    # now empty the rest of the buffer
    buf_size = len(qname_ind_map)
    while buf_size > 0:
        # buffer full so return first read
        yield buf[start_ind]                
        # delete the read qname to decrease the buffer size by
        # one and allow parser to iterate until it is full again
        return_qname = buf[start_ind][0][0].qname
        del qname_ind_map[return_qname]
        # advance start index
        start_ind += 1
        if start_ind == maxlen:
            start_ind = 0
        buf_size = len(qname_ind_map)

def sam_stdin_to_bam(output_bam_file, multihits):
    samfh = pysam.Samfile("-", "r")
    bamfh = pysam.Samfile(output_bam_file, "wb", template=samfh)
    for pe_reads in fix_alignment_order(samfh):
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
    sam_stdin_to_bam(sys.argv[1], int(sys.argv[2]))
