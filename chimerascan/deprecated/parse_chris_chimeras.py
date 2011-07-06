'''
Created on Jul 6, 2011

@author: mkiyer
'''
import sys
import subprocess

chrisfile = sys.argv[1]
seqfiles = sys.argv[2:4]


for line in open(chrisfile):
    fields = line.strip().split('\t')
    name = fields[0]    
    readids = fields[12:]
    fastq_lines = ([], [])
    for readid in readids:
        readid = readid[:-2]
        print 'chimera:', name, "read:", readid
        for mate,seqfile in enumerate(seqfiles):
            args = ["grep", "-m", "2", "-A", "2", readid, seqfile]
            res = subprocess.Popen(args, stdout=subprocess.PIPE).communicate()[0]
            fastq_lines[mate].extend(res.split('\n')[:-2])
    fh1 = open("%s_1.fq" % (name), "w")
    for line in fastq_lines[0]:
        print >>fh1, line
    fh1.close()
    fh2 = open("%s_2.fq" % (name), "w")
    for line in fastq_lines[1]:
        print >>fh2, line
    fh2.close()
