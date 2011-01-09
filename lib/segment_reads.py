'''
Created on Jan 7, 2011

@author: mkiyer
'''
def parse_fastq(line_iter):
    with line_iter:
        while True:
            lines = [line_iter.next().rstrip() for x in xrange(4)]
            yield lines

def parse_qname(qname, pe_trim=2):
    total_trim = pe_trim + 5
    newqname = qname[:-total_trim]
    mate = int(qname[-4])
    seg = int(qname[-3])
    num_segs = int(qname[-1])
    return newqname, mate, seg, num_segs

def segment_and_join_reads(infh1, infh2, outfh, segments):
    num_segs = len(segments)
    try:
        fastq_iters = [parse_fastq(infh1), parse_fastq(infh2)]
        while True:            
            for mate,fastq_iter in enumerate(fastq_iters):
                lines = fastq_iter.next()
                for seg_num,seg in enumerate(segments):
                    start, end = seg
                    newrname = "%s_%d%d/%d" % (lines[0][1:], mate, seg_num, num_segs)
                    newlines = ["@%s" % newrname,
                                lines[1][start:end],
                                "+%s" % newrname,
                                lines[3][start:end]]
                    print >>outfh, '\n'.join(newlines)
    except StopIteration:
        pass

if __name__ == '__main__':
    from optparse import OptionParser
    import sys
    parser = OptionParser("usage: %prog [options] <starts> <ends> <fastq>")
    options, args = parser.parse_args()
    starts = map(int, args[0].split(','))
    ends = map(int, args[1].split(','))
    segments = zip(starts, ends)
    assert len(segments) <= 10
    mate1_fastq_file = args[2]
    mate2_fastq_file = args[3]
    segment_and_join_reads(open(mate1_fastq_file), open(mate2_fastq_file), sys.stdout, segments)