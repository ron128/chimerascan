'''
Created on Jan 18, 2011

@author: mkiyer
'''
import argparse
import sys

def parse_fastq(line_iter):
    with line_iter:
        while True:
            lines = [line_iter.next().rstrip() for x in xrange(4)]
            yield lines

def parse_fasta(line_iter):    
    line = line_iter.next()
    assert line.startswith('>')
    tag = line.rstrip()[1:]
    seq = ''
    for line in line_iter:
        if line.startswith('>'):
            yield tag, seq
            tag = line.rstrip()[1:]
            seq = ''
        else:
            seq += line.rstrip()
    yield tag, seq

def segment_reads_fasta(infh, outfh, start, end):
    for tag,seq in parse_fasta(infh):
        newtag = ">%s_%03d-%03d" % (tag, start, end)
        print >>outfh, "%s\n%s" % (newtag, seq[start:end])

def segment_reads_fastq(infh, outfh, start, end):
    for lines in parse_fastq(infh):
        newrname = "%s_%03d-%03d" % (lines[0][1:], start, end)
        newlines = ["@%s" % newrname,
                    lines[1][start:end],
                    "+%s" % newrname,
                    lines[3][start:end]]
        print >>outfh, '\n'.join(newlines)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("sequence_file")
    parser.add_argument("start", type=int)
    parser.add_argument("end", type=int)
    options = parser.parse_args()
    segment_reads_fasta(open(options.sequence_file), sys.stdout, options.start, options.end)