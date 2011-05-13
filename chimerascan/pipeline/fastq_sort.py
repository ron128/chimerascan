'''
Created on Apr 28, 2011

@author: mkiyer
'''
import os
import tempfile
import logging

from chimerascan.lib.fastq_to_bam import fastq_to_bam, bam_to_fastq
from chimerascan.lib.config import JOB_SUCCESS, JOB_ERROR
from chimerascan import pysam


def sort_fastq_files(fastq_files, sorted_fastq_files, quals, tmp_dir):
    # convert to bam
    logging.debug("Converting FASTQ files to BAM")
    fd,tmpbam = tempfile.mkstemp(suffix=".bam", prefix="tmp", dir=tmp_dir)
    os.close(fd)
    fastq_to_bam(fastq_files, quals, tmpbam)
    # sort bam
    logging.debug("Sorting BAM")
    fd,srtbam = tempfile.mkstemp(suffix=".srt.bam", prefix="tmp", dir=tmp_dir)
    os.close(fd)
    pysam.sort("-n", tmpbam, os.path.splitext(srtbam)[0])
    # convert back to fastq
    logging.debug("Converting BAM to FASTQ")
    bam_to_fastq(srtbam, sorted_fastq_files)
    os.remove(tmpbam)
    os.remove(srtbam)
    return JOB_SUCCESS

if __name__ == '__main__':
    import sys
    sort_fastq_files(sys.argv[1:3], ["a1.fq", "a2.fq"], "solexa", ".")
