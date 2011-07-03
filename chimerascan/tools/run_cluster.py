#!/usr/bin/env python
'''
Created on Feb 1, 2011

@author: mkiyer

chimerascan: chimeric transcript discovery using RNA-seq

Copyright (C) 2011 Matthew Iyer

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''
import os
import subprocess

# TODO: change to allow customization
NODE_MEMORY = 45000.0
NODE_PROCESSORS = 12
MEM_PER_PROCESSOR = int(float(NODE_MEMORY) / NODE_PROCESSORS)
QUEUE_NAME = "flux"
ACCOUNT_NAME = "arul_flux"
HISEQ_WALLTIME = "80:00:00"
LOWSEQ_WALLTIME = "25:00:00"


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
             "#PBS -l qos=%s" % (ACCOUNT_NAME),
             "#PBS -A %s" % (ACCOUNT_NAME),
             "#PBS -q %s" % (QUEUE_NAME),
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
    from chimerascan import chimerascan_run    
    parser = chimerascan_run.RunConfig.get_option_parser()
    parser.set_usage("%prog [--big] <JOB_NAME> [chimerascan arguments]")
    parser.add_option("--big", dest="big", action="store_true", 
                      default=False, 
                      help="set this flag if you have a very large dataset "
                      "(more than 30M sequences) to adjust memory and "
                      "walltime limits accordingly")
    options, args = parser.parse_args()
    if options.big:
        walltime = HISEQ_WALLTIME
    else:
        walltime = LOWSEQ_WALLTIME
    job_name = args[0]    
    parser.remove_option("--big")
    if sys.argv.count("--big") > 0:
        sys.argv.remove("--big")    
    chimerascan_args = sys.argv[2:]
    # parse run parameters in config file and command line
    runconfig = chimerascan_run.RunConfig()
    runconfig.from_args(chimerascan_args, parser=parser)
    # now run chimerascan via 'qsub'
    args = [sys.executable,
            os.path.join(os.path.dirname(__file__),
                         chimerascan_run.__file__)] 
    args.extend(chimerascan_args)
    qsub(job_name, args, 
         runconfig.num_processors, 
         cwd=runconfig.output_dir, 
         walltime=walltime)
    sys.exit(0)

if __name__ == '__main__':
    main()
