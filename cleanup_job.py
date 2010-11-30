'''
Created on Nov 30, 2010

@author: mkiyer
'''
import sys
import argparse
import os
import logging
import subprocess
import shutil

from config import PipelineConfig, JobConfig, JOB_SUCCESS, JOB_ERROR

def cleanup_job(job_file, config_file):
    config = PipelineConfig.from_xml(config_file)
    job  = JobConfig.from_xml(job_file, config.output_dir)
    # remove temporary files that do not need to be copied
    if not os.path.exists(job.output_dir):
        logging.error("%s: Directory %s does not exist" % (job.name, job.output_dir))
        return JOB_ERROR
    for mate in xrange(len(job.fastq_files)):
        # Remove sequence files
        fastq_file = job.fastq_files[mate]
        if os.path.exists(fastq_file):
            logging.info("%s: Removing sequence file %s" % (job.name, fastq_file)) 
            os.remove(fastq_file)
    # see if this job is running on a remote system and needs to
    # retrieve/send data from/to a different server
    if job.remote:
        dst_output_dir = ":".join([job.remote_ip, job.dst_output_dir])
        copy_args = ["scp", "-r"]
    else:
        dst_output_dir = job.dst_output_dir
        copy_args = ["cp", "-r"]
    if job.output_dir == job.dst_output_dir:
        logging.info("%s: job output files are already in destination directory, so no need to copy results" % (job.name))
        return JOB_SUCCESS
    # copy the output data to the destination directory
    copy_args += [job.output_dir, dst_output_dir]
    logging.info("%s: Copying result directory %s to destination %s" % (job.name, job.output_dir, dst_output_dir))
    if subprocess.call(copy_args) != 0:
        logging.error("%s: Error copying result directory '%s'" % (job.name, str(copy_args)))
        return JOB_ERROR
    logging.info("%s: Copy successful" % (job.name))
    logging.info("%s: Deleting local files" % (job.name))
    #shutil.rmtree(job.output_dir)
    return JOB_SUCCESS

if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("job_file")
    parser.add_argument("config_file")
    options = parser.parse_args()
    sys.exit(cleanup_job(options.job_file, options.config_file))


