'''
Created on Oct 26, 2010

@author: mkiyer
'''
#Contains methods to tranform sequence strings
import string
#import lxml.etree as etree

def get_read_length(fastq_file):
    f = open(fastq_file)
    f.next()
    seq = f.next().strip()
    f.close()
    return len(seq)

def get_read_length_compressed(input_file):
    import gzip
    import bz2    
    import os
    suffix = os.path.splitext(input_file)[-1]
    if suffix == '.gz':
        f = gzip.GzipFile(input_file, 'r')
    elif suffix == '.bz2':
        f = bz2.BZ2File(input_file, 'r')
    else:
        f = open(input_file, 'r')
    f.next()
    seq = f.next().strip()
    f.close()
    return len(seq)

def parse_multihit_sam_file(samfh):    
    reads = []
    for read in samfh:        
        if len(reads) > 0 and read.qname != reads[-1].qname:
            yield reads
            reads = []
        reads.append(read)
    if len(reads) > 0:
        yield reads
        
CIGAR_M = 0 #match  Alignment match (can be a sequence match or mismatch)
CIGAR_I = 1 #insertion  Insertion to the reference
CIGAR_D = 2 #deletion  Deletion from the reference
CIGAR_N = 3 #skip  Skipped region from the reference
CIGAR_S = 4 #softclip  Soft clip on the read (clipped sequence present in <seq>)
CIGAR_H = 5 #hardclip  Hard clip on the read (clipped sequence NOT present in <seq>)
CIGAR_P = 6 #padding  Padding (silent deletion from the padded reference sequence)

def get_aligned_read_intervals(read):
    intervals = []
    # insert read into cluster tree
    astart,aend = read.pos, read.pos
    for op,length in read.cigar:
        if length == 0: continue
        if (op == CIGAR_I) or (op == CIGAR_S) or (op == CIGAR_H): continue
        if (op == CIGAR_P): assert False 
        if (op == CIGAR_N):
            assert astart != aend
            intervals.append((astart, aend))
            #print read.qname, read.cigar, ref, astart, aend
            astart = aend + length
        aend += length
    assert astart != aend
    if aend > astart:
        #print read.qname, read.cigar, ref, astart, aend
        intervals.append((astart, aend))
    assert aend == read.aend
    return intervals

def get_refs_from_bowtie_index(bowtie_index, split=True):
    import subprocess
    args = ['bowtie-inspect', '-s', bowtie_index]    
    p = subprocess.Popen(args, stdout=subprocess.PIPE)
    output = p.communicate()[0]
    refs = []
    for line in output.split('\n'):
        if not line:
            continue
        fields = line.split('\t')
        if fields[0].startswith('Sequence'):
            refname = fields[1].split()[0] if split else fields[1]
            refs.append((refname, int(fields[2])))
    return refs

class BEDGene():
    pass

def parse_bed12_line(line):
    if line is None:
        return None
    line = line.strip()
    if line.startswith('#'):
        return None
    if line.startswith('track'):
        return None
    fields = line.split('\t')
    # first six fields are required
    g = BEDGene()
    g.chrom = fields[0]
    g.tx_start = int(fields[1])
    g.tx_end = int(fields[2])
    g.name = fields[3]
    g.score = fields[4]
    g.strand = fields[5]        
    g.cds_start = int(fields[6])
    g.cds_end = int(fields[7])
    g.exon_count = int(fields[9])
    block_sizes = map(int, fields[10].split(',')[:-1])
    block_starts = map(int, fields[11].split(',')[:-1])        
    g.exon_starts = [(g.tx_start + start) for start in block_starts]        
    g.exon_ends = [(start + size) for start, size in zip(g.exon_starts, block_sizes)]
    g.exons = zip(g.exon_starts, g.exon_ends)
    g.introns = zip(g.exon_ends, g.exon_starts[1:])        
    return g

def parse_bed12_file(line_iter):
    '''
    parse a gene bed file
    '''
    for line in line_iter:
        g = parse_bed12_line(line)
        if g is None:
            continue
        yield g

#Translation table for reverse Complement, with ambiguity codes
DNA_COMPLEMENT = string.maketrans( "ACGTRYKMBDHVacgtrykmbdhv", "TGCAYRMKVHDBtgcayrmkvhdb" )
RNA_COMPLEMENT = string.maketrans( "ACGURYKMBDHVacgurykmbdhv", "UGCAYRMKVHDBugcayrmkvhdb" )
#Translation table for DNA <--> RNA
DNA_TO_RNA = string.maketrans( "Tt", "Uu" )
RNA_TO_DNA = string.maketrans( "Uu", "Tt" )

#reverse sequence string
def reverse( sequence ):
    return sequence[::-1]
#complement DNA sequence string
def DNA_complement( sequence ):
    return sequence.translate( DNA_COMPLEMENT )
#complement RNA sequence string
def RNA_complement( sequence ):
    return sequence.translate( RNA_COMPLEMENT )
#returns the reverse complement of the sequence
def DNA_reverse_complement( sequence ):
    sequence = reverse( sequence )
    return DNA_complement( sequence )
def RNA_reverse_complement( self, sequence ):
    sequence = reverse( sequence )
    return RNA_complement( sequence )
def to_DNA( sequence ):
    return sequence.translate( DNA_TO_RNA )
def to_RNA( sequence ):
    return sequence.translate( RNA_TO_DNA )



