'''
Created on Feb 1, 2011

@author: mkiyer
'''
import os
import subprocess

NODE_MEMORY = 45000.0
NODE_PROCESSORS = 12
MEM_PER_PROCESSOR = int(float(NODE_MEMORY) / NODE_PROCESSORS)
QUEUE_NAME = "arul_flux"

def qsub(job_name, args, num_processors, cwd=None, walltime="60:00:00", 
         pmem=None, deps=None, stdout=None, email_addresses=None):
    '''
    job_name: string name of job
    cmd: string (or list of args) containing command-line arguments for job
    num_processors: number of processors to submit with (cannot be greater than the number of processors per node)
    cwd: the "working directory" of the job (allows scripts to access files using relative pathnames)
    walltime: the walltime passed to qsub
    pmem: amount of memory allocated to this job
    deps: 'None' if no dependencies, of a python list of job ids
    stdout: string filename for storing stdout
    email_addresses: list of email addresses to send job information
    '''
    if isinstance(deps, basestring):
        deps = [deps]    
    if isinstance(email_addresses, basestring):
        email_addresses = [email_addresses]
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
             "#PBS -l qos=%s" % (QUEUE_NAME),
             "#PBS -A %s" % (QUEUE_NAME),
             "#PBS -q flux",
             "#PBS -V",
             "#PBS -j oe",
             "#PBS -o %s/%s" % (cwd, stdout)]
    if email_addresses is not None:        
        lines.append("#PBS -m bae")
        for email_address in email_addresses:
            lines.append("#PBS -M %s" % email_address)
    if deps is not None:
        lines.append("#PBS -W depend=afterok:%s" % (":".join([d for d in deps])))    
    
    if isinstance(args, basestring):
        lines.extend(["cat $PBS_NODEFILE", 
                  "cd %s" % (cwd),
                  args])
    else:
        lines.extend(["cat $PBS_NODEFILE", 
                  "cd %s" % (cwd), 
                  ' '.join(args)])
    print '\n'.join(lines)
    p = subprocess.Popen("qsub", stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    p.stdin.write('\n'.join(lines))
    job_id = p.communicate()[0]
    #job_id='test'
    return job_id.strip()

def main():
    import sys
    from run_chimerascan import RunConfig
    job_name = sys.argv[1]
    # parse run parameters in config file and command line
    runconfig = RunConfig()
    runconfig.from_args(sys.argv[2:])
    args = [sys.executable,
            os.path.join(os.path.dirname(__file__), 
                         "run_chimerascan.py")]
    args.extend(sys.argv[1:])
    qsub(job_name, args, 
         runconfig.num_processors, 
         cwd=runconfig.output_dir, 
         walltime="60:00:00")
    sys.exit(0)

if __name__ == '__main__':
    main()
