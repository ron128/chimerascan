'''
Created on Jan 7, 2011

@author: mkiyer
'''
import re

# currently we add 5 characters to the end of each rname to
# encode the mate, segment, and number of segments in the 
# read
SUFFIX_LENGTH = 5

def parse_fastq(line_iter):
    with line_iter:
        while True:
            lines = [line_iter.next().rstrip() for x in xrange(4)]
            yield lines

def parse_qname(qname):
    newqname = qname[:-SUFFIX_LENGTH]
    mate = int(qname[-4])
    seg = int(qname[-3])
    num_segs = int(qname[-1])
    return newqname, mate, seg, num_segs

def segment_and_join_reads(infhs, outfh, segments):
    num_segs = len(segments)
    # make regular expression to clip the trailing '/1' and '/2'
    # from paired end reads
    suffix_length = None
    try:
        fastq_iters = [parse_fastq(fh) for fh in infhs]
        # now suffix lengths are known and can automatically remove them 
        # without checking first
        while True:            
            for mate,fastq_iter in enumerate(fastq_iters):
                lines = fastq_iter.next()
                qname = lines[0][1:]
                if suffix_length is None:
                    # parse first read and find suffix length
                    newqname = re.split(r'/\d$', qname)[0]
                    suffix_length = len(qname) - len(newqname)                    
                for seg_num,seg in enumerate(segments):
                    start, end = seg
                    newqname = "%s_%d%d/%d" % (qname[:len(qname)-suffix_length], mate, seg_num, num_segs)
                    newlines = ["@%s" % newqname,
                                lines[1][start:end],
                                "+%s" % newqname,
                                lines[3][start:end]]
                    print >>outfh, '\n'.join(newlines)
    except StopIteration:
        pass

if __name__ == '__main__':
    from optparse import OptionParser
    import sys
    parser = OptionParser("usage: %prog [options] <starts> <ends> <fastq_file> [<fastq_file>]")
    options, args = parser.parse_args()
    starts = map(int, args[0].split(','))
    ends = map(int, args[1].split(','))
    segments = zip(starts, ends)
    assert len(segments) <= 10
    fastq_files = args[2:]
    fastq_fhs = [open(f) for f in fastq_files]
    segment_and_join_reads(fastq_fhs, sys.stdout, segments)