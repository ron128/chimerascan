'''
Created on Apr 28, 2011

@author: mkiyer
'''
from chimerascan import pysam
from math import log10
from string import maketrans

def get_solexa_qual_conversion_table():
    """
    return a translation table that can be used by str.translate() for
    converting solexa to sanger quality scores
    """
    offset = 64
    conv_table = ['!'] * 256
    conv_table[offset:] = "I" * (256-offset)
    for solq in xrange(-5, 40):
        phredq = 10*log10(1 + 10**(solq/10.0))
        phredchr = chr(int(round(33 + phredq)))
        conv_table[offset + solq] = phredchr
    conv_string = ''.join(conv_table)
    return maketrans(''.join(map(chr, range(256))), conv_string)

def get_illumina_qual_conversion_table():
    """Illumina 1.3+ format can encode a Phred quality score from 0 to 62 
    using ASCII 64 to 126 (although in raw read data Phred scores from 0 
    to 40 only are expected).
    """
    offset = 64
    conv_table = ['!'] * 256
    for x in xrange(0, 62):
        conv_table[offset+x] = chr(33 + x)
    conv_table[offset+40:] = "I" * (256-(offset+40))
    conv_string = ''.join(conv_table)
    return maketrans(''.join(map(chr, range(256))), conv_string)    

def get_sanger_qual_conversion_table():
    offset = 33
    tbl = map(chr, range(256))
    tbl[:offset] = "!" * offset
    tbl[offset+40:] = "I" * (256-(offset+40))
    return maketrans(''.join(map(chr, range(256))), ''.join(tbl))

conv_tables = {"sanger": get_sanger_qual_conversion_table(),
               "illumina": get_illumina_qual_conversion_table(),
               "solexa": get_solexa_qual_conversion_table()}

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
    qual_trans_table = conv_tables[qual_format]
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
                a.qual = qual.translate(qual_trans_table)
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
    sol2std = get_solexa_qual_conversion_table()
    illumina2std = get_illumina_qual_conversion_table()
    import sys
    fastq_to_bam(["read1.fq", "read2.fq"], "solexa", "hi.bam")
    bam_to_fastq("hi.bam", ["a1.fq", "a2.fq"])
