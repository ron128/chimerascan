'''
Created on Apr 28, 2011

@author: mkiyer
'''
from chimerascan import pysam
from seq import get_qual_conversion_func

def parse_fastq(line_iter):
    with line_iter:
        while True:
            rid = line_iter.next().rstrip()[1:]
            seq = line_iter.next().rstrip()
            line_iter.next()
            qual = line_iter.next().rstrip()
            yield rid, seq, qual

def fastq_to_bam(fastq_files, qual_format, bam_file):
    fqfhs = [parse_fastq(open(f)) for f in fastq_files]
    qual_func = get_qual_conversion_func(qual_format)
    header = {'HD': {'VN': '1.0', 'SO': 'unknown'}}
#              'SQ': [{'LN': 1, 'SN': 'dummy'}]}
    bamfh = pysam.Samfile(bam_file, "wb", header=header)    
    try:
        while True:
            for i,fqiter in enumerate(fqfhs):
                id,seq,qual = fqiter.next()
                a = pysam.AlignedRead()
                a.rname = -1
                a.mrnm = -1
                #a.pos = 0
                #a.mpos = 0
                a.qname = id
                a.seq = seq
                a.qual = qual_func(qual)
                a.is_read1 = (i == 0)
                a.is_read2 = (i == 1)
                bamfh.write(a)
    except StopIteration:
        pass
    bamfh.close()  

def bam_to_fastq(bam_file, fastq_files):
    fqfhs = [open(f, "w") for f in fastq_files]
    bamfh = pysam.Samfile(bam_file, "rb")
    for r in bamfh:
        if r.is_read1:
            i = 0
        elif r.is_read2:
            i = 1
        record = "@%s\n%s\n+\n%s" % (r.qname,r.seq,r.qual)
        print >>fqfhs[i], record

if __name__ == '__main__':
    sol2std = get_qual_conversion_func("solexa")
    illumina2std = get_qual_conversion_func("illumina")
    import sys
    fastq_to_bam(["read1.fq", "read2.fq"], "solexa", "hi.bam")
    bam_to_fastq("hi.bam", ["a1.fq", "a2.fq"])
