'''
Created on Oct 26, 2010

@author: mkiyer
'''

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

def parse_bed12_line(line):
    class BEDGene():
        pass
    if line is None:
        return None
    line = line.strip()
    if line.startswith('#'):
        logging.debug("skipping comment line: %s" % (line))
        return None
    if line.startswith('track'):
        logging.debug("skipping track header line: %s"  % (line))
        return None
    thisfields = line.split('\t')
    # first six fields are required
    g = BEDGene()
    g.chrom = thisfields[0]
    g.tx_start = int(thisfields[1])
    g.tx_end = int(thisfields[2])
    g.name = thisfields[3]
    if len(thisfields) <= 4:
        g.score = 0
        g.strand = '.'
    else:
        g.score = thisfields[4]
        g.strand = thisfields[5]        
    if len(thisfields) <= 6:
        g.cds_start = g.tx_start
        g.cds_end = g.tx_end
        g.exon_count = 1
        g.exons = [(g.tx_start, g.tx_end)]
        g.introns = []
    else:
        g.cds_start = int(thisfields[6])
        g.cds_end = int(thisfields[7])
        g.exon_count = int(thisfields[9])
        block_sizes = map(int, thisfields[10].split(',')[:-1])
        block_starts = map(int, thisfields[11].split(',')[:-1])        
        g.exon_starts = [(g.tx_start + start) for start in block_starts]        
        g.exon_ends = [(start + size) for start, size in zip(g.exon_starts, block_sizes)]
        g.exons = zip(g.exon_starts, g.exon_ends)
        g.introns = zip(g.exon_ends, g.exon_starts[1:])        
    return g


#Contains methods to tranform sequence strings
import string

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
