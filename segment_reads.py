'''
Created on Oct 14, 2010

@author: mkiyer
'''
import argparse
import sys

def parse_fastq(line_iter):
    with line_iter:
        while True:
            lines = [line_iter.next().rstrip() for x in xrange(4)]
            yield lines

def segment_reads(infh, outfh, start, end):
    for lines in parse_fastq(infh):
        newrname = "%s_%03d-%03d" % (lines[0][1:], start, end)
        newlines = ["@%s" % newrname,
                    lines[1][start:end],
                    "+%s" % newrname,
                    lines[3][start:end]]
        print >>outfh, '\n'.join(newlines)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("fastq_file")
    parser.add_argument("start", type=int)
    parser.add_argument("end", type=int)
    options = parser.parse_args()
    segment_reads(open(options.fastq_file), sys.stdout, options.start, options.end)