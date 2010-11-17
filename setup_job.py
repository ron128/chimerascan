'''
Created on Nov 17, 2010

@author: mkiyer
'''
import sys
import os
import argparse
import logging
import subprocess

from config import PipelineConfig, JobConfig, JOB_SUCCESS, JOB_ERROR

def copy_file(input_file, output_file, remote=False, remote_ip=None):
    if remote:
        input_file = ":".join(remote_ip, input_file)
        args = ["scp", input_file, output_file]
    else:
        args = ["cp", input_file, output_file]
    if subprocess.call(args) != 0:
        logging.error("Error copying sequence file args='%s'" % str(args))
        return 1
    return 0

def decompress_file(input_file, output_file):
    suffix = os.path.splitext(input_file)[-1]
    if suffix == '.gz':
        args = ["zcat", input_file]
    elif suffix == '.bz2':
        args = ["bzcat", input_file]
    else:
        args = ["cat", input_file]
    if subprocess.call(args, stdout=open(output_file, "wb")) != 0:
        logging.error("Error decompressing file args='%s'" % str(args))
        return 1
    return 0

def copy_sequence_job(job_file, config_file):
    config = PipelineConfig.from_xml(config_file)    
    job  = JobConfig.from_xml(job_file, config.output_dir)
    # make output directory
    if not os.path.exists(job.output_dir):
        os.makedirs(job.output_dir)
        logging.info("%s: created output directory %s" % (job.name, job.output_dir))
    # copy fastq files    
    for mate in xrange(len(job.src_fastq_files)):
        src_fastq_file = job.src_fastq_files[mate]
        dst_fastq_file = job.dst_fastq_files[mate]
        fastq_file = job.fastq_files[mate]
        # do not copy if file exists
        if os.path.exists(fastq_file):
            logging.info("Local fastq file %s exists, skipping copy" % dst_fastq_file)
            continue
        # copy sequence files to local output dir
        logging.info("%s: Copying sequence file %s to %s remote=%s" % 
                     (job.name, src_fastq_file, dst_fastq_file, str(job.remote)))
        if copy_file(src_fastq_file, dst_fastq_file, job.remote, job.remote_ip) != 0:
            logging.error("Error copying file %s" % dst_fastq_file)
            return JOB_ERROR
        # uncompress if necessary
        logging.info("%s: Uncompressing file %s to %s" %
                     (job.name, dst_fastq_file, fastq_file)) 
        if decompress_file(dst_fastq_file, fastq_file) != 0:
            logging.error("Error decompressing file %s" % dst_fastq_file)
            return JOB_ERROR
        # remote tmp fastq file intermediate
        os.remove(dst_fastq_file)
    return JOB_SUCCESS

if __name__ == '__main__': 
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("config_file")
    parser.add_argument("job_file")
    options = parser.parse_args()
    sys.exit(copy_sequence_job(options.job_file, options.config_file))
