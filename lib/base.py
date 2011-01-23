'''
Created on Oct 26, 2010

@author: mkiyer
'''
import os
import string
import subprocess

def parse_library_type(library_type):
    s1 = 0 if library_type[0] == 'f' else 1
    s2 = 0 if library_type[1] == 'f' else 1
    return (s1, s2)

def get_read_length(fastq_file):
    f = open(fastq_file)
    f.next()
    seq = f.next().strip()
    f.close()
    return len(seq)

def get_read_length_compressed(input_file):
    import gzip
    import bz2    
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

def check_executable(filename):
    # check that samtools binary exists
    devnullfh = open(os.devnull, 'w')        
    try:
        subprocess.call([filename], stdout=devnullfh, stderr=devnullfh)
    except OSError:
        return False
    devnullfh.close()
    return True

def parse_multihit_sam_file(samfh):    
    reads = []
    for read in samfh:        
        if len(reads) > 0 and read.qname != reads[-1].qname:
            yield reads
            reads = []
        reads.append(read)
    if len(reads) > 0:
        yield reads

def parse_multihit_alignments(samfh):
    buf = []
    ind = 0
    for read in samfh:
        if (ind > 0) and (read.qname != buf[ind-1].qname):
            yield buf[:ind]
            ind = 0
        if ind < len(buf):
            buf[ind] = read
        else:
            buf.append(read)
        ind += 1
    if ind > 0:
        yield buf[:ind]

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




