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

from config import PipelineConfig, JobConfig, JOB_SUCCESS, JOB_ERROR

_module_dir = os.path.abspath(os.path.dirname(__file__))

def run_chimerascan_pipeline(task):
    job_file, config_file = task
    logging.info("Job %s starting.." % (job_file))    
    config = PipelineConfig.from_xml(config_file)    
    job = JobConfig.from_xml(job_file, config.output_dir)
    # setup job
    if not os.path.exists(job.output_dir):
        os.makedirs(job.output_dir)
        logging.info("%s: created output directory %s" % (job.name, job.output_dir))    
    logging.info("%s: Setting up job at output dir %s" % (job.name, job.output_dir))    
    py_script = os.path.join(_module_dir, "setup_job.py")
    args = [sys.executable, py_script, config_file, job_file]
    fout = open(os.path.join(job.output_dir, "setup.log"), "w")
    retcode = subprocess.call(args, stderr=fout)
    fout.close()
    if retcode != 0:
        logging.error("%s: error setting up job" % (job.name))    
        return job_file, JOB_ERROR
    # make output directory
    if not os.path.exists(job.chimerascan_dir):
        os.makedirs(job.chimerascan_dir)
        logging.info("%s: created output directory %s" % (job.name, job.chimerascan_dir))
    # align reads and find discordant mappings
    logging.info("%s: Aligning reads" % (job.name))    
    py_script = os.path.join(_module_dir, "align.py")
    args = [sys.executable, py_script,
            "--bowtie-bin", config.bowtie_bin,
            "--bowtie-index", config.bowtie_index,
            "--bowtie-threads", config.bowtie_threads,
            "--multihits", config.multihits,
            "--mismatches", config.mismatches,
            "--seed-length", config.seed_length,
            "--gene-bed", config.gene_bed_file,
            "--gene-fasta-prefix", config.gene_fasta_prefix,
            "--quals", job.fastq_format,
            job.fastq_files[0],
            job.fastq_files[1],
            job.discordant_bam_file,
            job.expression_file]
    args = map(str, args)
    logging.debug("Running alignment with args: %s" % ' '.join(args))
    fout = open(os.path.join(job.output_dir, "align.log"), "w")
    retcode = subprocess.call(args, stderr=fout)
    fout.close()
    if retcode != 0:
        logging.error("%s: error aligning reads" % (job.name))    
        return job_file, JOB_ERROR    
    return job_file, JOB_SUCCESS
    # nominate chimeras
    
    # convert chimeras to fasta file
    # build bowtie index from fasta file
    # map reads to the new index
    # integrate spanning and encompassing reads


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
    imap_unordered_it = pool.imap_unordered(run_chimerascan_pipeline, tasks)
    for res in imap_unordered_it:
        job_file, retcode = res
        logging.info("Job %s finished with code %d" % (job_file, retcode))
    pool.close()
    pool.join()

if __name__ == '__main__': main()