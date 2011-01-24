'''
Created on Nov 17, 2010

@author: mkiyer
'''
import sys
import os
import logging
import argparse
import subprocess
import lxml.etree as etree

from config import PipelineConfig, JobConfig, JOB_SUCCESS, JOB_ERROR
from setup_job import copy_sequence_job
from base import get_read_length_compressed

_module_dir = os.path.abspath(os.path.dirname(__file__))

NODE_MEMORY = 45000.0
NODE_PROCESSORS = 12
MEM_PER_PROCESSOR = int(float(NODE_MEMORY) / NODE_PROCESSORS)

def qsub(job_name, cmd, num_processors, cwd=None, walltime="60:00:00", pmem=None, deps=None, stdout=None, email=False):
    if cwd is None:
        cwd = os.getcwd()
    num_processors = min(NODE_PROCESSORS, num_processors)
    if pmem is None:
        pmem = MEM_PER_PROCESSOR
    if stdout is None:
        stdout = "${PBS_JOBID}"
    lines = ["#!/bin/sh",
             "#PBS -N %s" % job_name,
             "#PBS -l nodes=%d:ppn=%d,walltime=%s,pmem=%dmb" % (1, num_processors, walltime, pmem),
             "#PBS -l qos=arul_flux",
             "#PBS -A arul_flux",
             "#PBS -q flux",
             "#PBS -V",
             "#PBS -j oe",
             "#PBS -o %s/%s" % (cwd, stdout)]
    if email:
        lines.extend(["#PBS -m ae",
                      "#PBS -M chrmaher@umich.edu",
                      "#PBS -M mkiyer@umich.edu",
                      "#PBS -M trbarret@med.umich.edu"])        
    if deps is not None:
        lines.append("#PBS -W depend=afterok:%s" % (":".join([d for d in deps])))    
    lines.extend(["cat $PBS_NODEFILE", 
                  "cd %s" % (cwd), 
                  cmd])
    print '\n'.join(lines)
    p = subprocess.Popen("qsub", stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    p.stdin.write('\n'.join(lines))
    job_id = p.communicate()[0]
    return job_id.strip()

def up_to_date(outfile, infile):
    if not os.path.exists(infile):
        return False
    if not os.path.exists(outfile):
        return False
    if os.path.getsize(outfile) == 0:
        return False    
    return os.path.getmtime(outfile) >= os.path.getmtime(infile)

def run_job_on_cluster(job_file, config_file):
    logging.info("Job %s starting.." % (job_file))    
    config = PipelineConfig.from_xml(config_file)    
    job = JobConfig.from_xml(job_file, config.output_dir)
    deps = []
    #
    # Setup job by copying sequences and uncompressing them
    #
    if not os.path.exists(job.output_dir):
        os.makedirs(job.output_dir)
        logging.info("%s: Created output directory %s" % (job.name, job.output_dir))    
    logging.info("%s: Setting up job at output dir %s" % (job.name, job.output_dir))    
    retcode = copy_sequence_job(job_file, config_file)
    if retcode != 0:
        logging.error("%s: Error setting up job" % (job.name))    
        return JOB_ERROR
    #
    # Get the read length from the fastq files
    #
    if os.path.exists(job.dst_fastq_files[0]):
        read_length = get_read_length_compressed(job.dst_fastq_files[0])
    elif os.path.exists(job.fastq_files[0]):
        read_length = get_read_length_compressed(job.fastq_files[0])
    logging.info("%s: First sequence in file has read length %d" % (job.name, read_length))    
    #
    # Uncompress sequences
    #
    if all(up_to_date(job.fastq_files[mate], job.src_fastq_files[mate]) for mate in xrange(len(job.fastq_files))):
        logging.info("[SKIPPED] Uncompressed sequence files %s are up to date" % (job.fastq_files))
    else:
        logging.info("%s: Uncompressing sequences" % (job.name))    
        py_script = os.path.join(_module_dir, "setup_job.py")
        args = [sys.executable, py_script, "--uncompress", config_file, job_file]
        cmd = ' '.join(args)
        job_id = qsub(job.name, cmd, 1, cwd=job.output_dir, walltime="10:00:00", stdout="uncompress.log", email=False)
        deps = [job_id]
    #
    # Discordant reads alignment 
    #
    # make output directory
    if not os.path.exists(job.chimerascan_dir):
        os.makedirs(job.chimerascan_dir)
        logging.info("%s: created output directory %s" % (job.name, job.chimerascan_dir))
    if all(up_to_date(job.discordant_bam_file, f) for f in job.fastq_files):
        logging.info("[SKIPPED] Discordant reads alignment %s is up to date" % (job.discordant_bam_file))
    else:
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
        cmd = ' '.join(args)
        num_processors = NODE_PROCESSORS
        job_id = qsub(job.name, cmd, num_processors, cwd=job.output_dir, walltime="40:00:00", deps=deps, stdout="align.log", email=False)
        deps = [job_id]
    #
    # Nominate chimeras
    #
    if up_to_date(job.chimera_bedpe_file, job.discordant_bam_file):
        logging.info("[SKIPPED] Nominate chimeras file %s is up to date" % (job.chimera_bedpe_file))
    else:
        py_script = os.path.join(_module_dir, "nominate_chimeras.py")
        args = [sys.executable, py_script,
                "--bedtools-path", config.bedtools_path,
                "--gene-bed", config.gene_bed_file,
                "--gene-name", config.gene_name_file,
                job.name,
                job.discordant_bam_file,
                job.chimerascan_dir,
                job.chimera_bedpe_file]
        args = map(str, args)
        cmd = ' '.join(args)
        job_id = qsub(job.name, cmd, num_processors=1, cwd=job.output_dir, walltime="40:00:00", deps=deps, stdout="chimeras.log", email=False)
        deps = [job_id]
    #
    # Convert BEDPE to FASTA format
    #
    if up_to_date(job.chimera_fasta_file, job.chimera_bedpe_file):
        logging.info("[SKIPPED] BEDPE to FASTA conversion file %s is up to date" % (job.chimera_fasta_file))
    else:
        py_script = os.path.join(_module_dir, "bedpe_to_fasta.py")
        args = [sys.executable, py_script,
                "--rlen", read_length,
                "--gene-fasta-prefix", config.gene_fasta_prefix,
                job.chimera_bedpe_file,
                config.ref_fasta_file,
                job.chimera_fasta_file,
                job.chimera_mapping_file]              
        args = map(str, args)
        cmd = ' '.join(args)
        job_id = qsub(job.name, cmd, num_processors=1, cwd=job.output_dir, walltime="1:00:00", deps=deps, stdout="bedpe_to_fasta.log", email=False)
        deps = [job_id]
    # 
    # Build bowtie index
    #
    bowtie_chimera_index_file = job.bowtie_chimera_index + ".1.ebwt"
    if up_to_date(bowtie_chimera_index_file, job.chimera_fasta_file):
        logging.info("[SKIPPED] Bowtie index %s is up to date" % (job.bowtie_chimera_index))
    else:
        args = [config.bowtie_build, job.chimera_fasta_file, job.bowtie_chimera_index]
        cmd = ' '.join(args)
        job_id = qsub(job.name, cmd, num_processors=1, cwd=job.output_dir, walltime="1:00:00", deps=deps, stdout="bowtie_build.log", email=False)
        deps = [job_id]
    #
    # Spanning read alignment
    #
    py_script = os.path.join(_module_dir, "segmented_spanning_align.py")
    job_ids = []
    for mate,fq in enumerate(job.fastq_files):
        if up_to_date(job.spanning_bowtie_output_files[mate], bowtie_chimera_index_file):
            logging.info("[SKIPPED] Spanning read alignment file %s is up to date" % (job.spanning_bowtie_output_files[mate]))
            continue
        args = [sys.executable, py_script,
                "--bowtie-bin", config.bowtie_bin,
                "--bowtie-index", job.bowtie_chimera_index,
                "--bowtie-threads", 2,
                "--segment-multihits", 2,
                "--segment-mismatches", 2,
                "--segment-length", 25,
                "--quals", job.fastq_format,
                fq,
                job.spanning_bowtie_output_files[mate]]
        cmd = ' '.join(map(str, args))
        job_ids.append(qsub(job.name, cmd, num_processors=NODE_PROCESSORS, cwd=job.output_dir, walltime="20:00:00", deps=deps, 
                            stdout="spanning_mate%d.log" % mate, email=False))
    #
    # Synthesis of spanning and encompassing reads
    #
    if all(up_to_date(job.spanning_chimera_file, f) for f in job.spanning_bowtie_output_files):
        logging.info("[SKIPPED] Processed spanning alignment %s is up to date" % (job.spanning_chimera_file))
    else:
        py_script = os.path.join(_module_dir, "process_spanning_alignments.py")
        args = [sys.executable, py_script,
                "--rlen", read_length,
                "--anchor-min", config.anchor_min,
                "--anchor-max", config.anchor_max,
                "--anchor-mismatches", config.anchor_mismatches,
                job.chimera_mapping_file,
                job.spanning_chimera_file] + job.spanning_bowtie_output_files
        cmd = ' '.join(map(str, args))
        qsub(job.name, cmd, num_processors=1, cwd=job.output_dir, walltime="2:00:00", deps=job_ids, 
             stdout="process_spanning_alignments.log", email=True)
    return JOB_SUCCESS

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    #parser.add_argument("-p", "--processes", type=int, dest="num_processors", default=1)
    parser.add_argument("job_file")
    parser.add_argument("config_file")
    options = parser.parse_args()
    run_job_on_cluster(options.job_file, options.config_file)
    
if __name__ == '__main__': main()