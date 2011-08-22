'''
Created on Jan 5, 2011

@author: Dan Blankenberg

Code from the Galaxy project (http://galaxy.psu.edu)
Contains methods to transform sequence strings
'''
import string
from math import log10
from string import maketrans

# Quality score formats
SANGER_FORMAT = "sanger"
SOLEXA_FORMAT = "solexa"
ILLUMINA_FORMAT = "illumina"
FASTQ_QUAL_FORMATS = [SANGER_FORMAT, SOLEXA_FORMAT, ILLUMINA_FORMAT]

#Translation table for reverse Complement, with ambiguity codes
DNA_COMPLEMENT = string.maketrans( "ACGTRYKMBDHVacgtrykmbdhv", "TGCAYRMKVHDBtgcayrmkvhdb" )
RNA_COMPLEMENT = string.maketrans( "ACGURYKMBDHVacgurykmbdhv", "UGCAYRMKVHDBugcayrmkvhdb" )
#Translation table for DNA <--> RNA
DNA_TO_RNA = string.maketrans( "Tt", "Uu" )
RNA_TO_DNA = string.maketrans( "Uu", "Tt" )

def DNA_complement( sequence ):
    '''complement DNA sequence string'''
    return sequence.translate( DNA_COMPLEMENT )
def DNA_reverse_complement( sequence ):
    '''returns the reverse complement of the sequence'''
    return DNA_complement(sequence[::-1])
def to_DNA( sequence ):
    return sequence.translate( DNA_TO_RNA )
#complement RNA sequence string
def RNA_complement( sequence ):
    return sequence.translate( RNA_COMPLEMENT )
def RNA_reverse_complement( self, sequence ):
    return RNA_complement( sequence[::-1] )
def to_RNA( sequence ):
    return sequence.translate( RNA_TO_DNA )

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

def get_qual_conversion_func(qual_format):
    conv_tables = {SANGER_FORMAT: get_sanger_qual_conversion_table(),
                   ILLUMINA_FORMAT: get_illumina_qual_conversion_table(),
                   SOLEXA_FORMAT: get_solexa_qual_conversion_table()}
    tbl = conv_tables[qual_format]
    return lambda q: q.translate(tbl)

class FASTQRecord:
    __slots__ = ("qname", "seq", "qual", "readnum")
    def __init__(self, qname, seq, qual, readnum):
        self.qname = qname
        self.seq = seq
        self.qual = qual
        self.readnum = readnum
        
    def to_string(self):
        return ("@%s/%d\n%s\n+\n%s" % 
                (self.qname, self.readnum, self.seq, self.qual))

def parse_fastq_record(line_iter, 
                       convert_quals=False,
                       qual_format=SANGER_FORMAT):
    qual_func = get_qual_conversion_func(qual_format)    
    try:        
        qname = line_iter.next().rstrip()[1:]
        readnum = int(qname[-1])
        qname = qname[:-2]
        seq = line_iter.next().rstrip()
        line_iter.next()
        qual = line_iter.next().rstrip()
        if convert_quals:
            qual = qual_func(qual)
        yield FASTQRecord(qname, seq, qual, readnum)
        while True:
            # qname
            qname = line_iter.next().rstrip()[1:]
            readnum = int(qname[-1])
            qname = qname[:-2]
            # seq
            seq = line_iter.next().rstrip()
            # qname again (skip)
            line_iter.next()
            # qual
            qual = line_iter.next().rstrip()
            if convert_quals:
                qual = qual_func(qual)
            yield FASTQRecord(qname, seq, qual, readnum)
    except StopIteration:
        pass

def calc_homology(seq1, seq2, num_mismatches):
    smallest_len = min(len(seq1), len(seq2))
    mm = 0
    i = 0
    for i in xrange(smallest_len):
        if seq1[i] != seq2[i]:
            mm += 1
            if mm > num_mismatches:
                return i
    return i + 1

BASES_PER_LINE = 50
def split_seq(seq, chars_per_line=BASES_PER_LINE):
    pos = 0
    newseq = []
    while pos < len(seq):
        if pos + chars_per_line > len(seq):        
            endpos = len(seq)
        else:
            endpos = pos + chars_per_line
        newseq.append(seq[pos:endpos])
        pos = endpos
    return '\n'.join(newseq)
