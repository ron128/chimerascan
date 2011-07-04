'''
Created on Jul 3, 2011

@author: mkiyer
'''

def parse_fastq(line_iter):
    with line_iter:
        while True:
            lines = [line_iter.next().rstrip() for x in xrange(4)]
            yield lines

def main():
    from optparse import OptionParser
    parser = OptionParser("usage: %prog [options] <in.fq>")
    options, args = parser.parse_args()
    infq = args[0]
    for lines in parse_fastq(open(infq)):
        qname_string = lines[0]
        qname,info = qname_string.split()
        readnum = int(info[0])
        lines[0] = "%s/%d" % (qname,readnum)
        lines[2] = "+"
        print '\n'.join(lines)

if __name__ == '__main__':
    main()
