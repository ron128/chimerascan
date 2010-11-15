'''
Created on Oct 26, 2010

@author: mkiyer
'''
import os
import sys
import argparse
import multiprocessing
import subprocess
import logging

from base import PipelineConfig, JobConfig, JOB_SUCCESS, JOB_ERROR

_module_dir = os.path.abspath(os.path.dirname(__file__))

def run_fusion_pipeline(task):
    # align reads and find discordant mappings
    
    # cluster and find fusion coverage islands
    
    # splice coverage islands together and build fasta file
    
    # build a bowtie index

    # map the initially unmapped reads to the new index
    
    # count fusion spanning reads and nominate fusions

    job_file, config_file = task
    logging.info("Job %s starting.." % (job_file))    
    config = PipelineConfig.from_xml(config_file)    
    job = JobConfig.from_xml(job_file, config.output_dir)
    if not os.path.exists(job.output_dir):
        logging.error("%s: missing output directory %s" % (job.name, job.output_dir))    
        return job_file, 1

    logging.info("%s: Creating UCSC tracks" % (job.name))    
    py_script = os.path.join(_module_dir, "make_ucsc_tracks.py")
    args = [sys.executable, py_script, job_file, config_file]
    fout = open(os.path.join(job.output_dir, "ucsc.out"), "w")
    retcode = subprocess.call(args, stderr=fout)
    fout.close()
    if retcode != 0:
        logging.error("%s: error creating UCSC tracks" % (job.name))    
        return job_file, 1
    
    logging.info("%s: Creating pileup tracks" % (job.name))    
    py_script = os.path.join(_module_dir, "make_pileup_track.py")
    args = [sys.executable, py_script, job_file, config_file]
    fout = open(os.path.join(job.output_dir, "pileup.out"), "w")
    retcode = subprocess.call(args, stderr=fout)
    fout.close()
    if retcode != 0:
        logging.error("%s: error creating pileup tracks" % (job.name))    
        return job_file, 1
    return job_file, 0

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--processes", type=int, dest="num_processors", default=1)
    parser.add_argument("config_file")
    parser.add_argument("job_files", nargs="+")
    options = parser.parse_args()
    tasks = []
    for job_file in options.job_files:
        tasks.append((job_file, options.config_file))
    # start worker processes
    pool = multiprocessing.Pool(options.num_processors)
    imap_unordered_it = pool.imap_unordered(run_job, tasks)
    for res in imap_unordered_it:
        job_file, retcode = res
        logging.info("Job %s finished with code %d" % (job_file, retcode))
    pool.close()
    pool.join()

if __name__ == '__main__': main()