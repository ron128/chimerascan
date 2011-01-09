'''
Created on Jan 8, 2011

@author: mkiyer
'''

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


def main():
    from optparse import OptionParser
    parser = OptionParser("usage: %prog [options] <seq.fa> <outprefix>")
    parser.add_option("--rlen", dest="read_length", type="int", default=50)
    parser.add_option("--isize", dest="insert_size", type="int", default=200)
    options, args = parser.parse_args()
    # extract command line arguments
    fasta_file = args[0]
    outprefix = args[1]
    fasta_records = list(parse_fasta(open(fasta_file)))
    tag, seq = fasta_records[0]
    tag = tag.split()[0]
    read_length = options.read_length
    insert_size = options.insert_size
    assert insert_size > read_length
    read1 = (0, read_length)
    read2 = (insert_size - read_length, insert_size)

    fh1 = open(outprefix + '1.fq', 'w')
    fh2 = open(outprefix + '2.fq', 'w')
    for end in xrange(insert_size, len(seq)):
        start = end - insert_size        
        start1,end1 = start + read1[0], start + read1[1]
        start2,end2 = start + read2[0], start + read2[1]        
        seq1 = seq[start1:end1]
        seq2 = seq[start2:end2]
        thistag = '%s:%d-%d' % (tag, start1, end2)
        print >>fh1, '@%s/1\n%s\n+%s/1\n%s' % (thistag, seq1, thistag, 'a'*read_length)        
        print >>fh2, '@%s/2\n%s\n+%s/2\n%s' % (thistag, seq2, thistag, 'a'*read_length)
    fh1.close()
    fh2.close()

#@PATHBIO-SOLEXA2_30TUEAAXX:1:3:1574:973/1
#ATAAGTACAGTCTATTGATCCGTGATGCATGCTGTACATGAATTTGGTGNNNN
#+PATHBIO-SOLEXA2_30TUEAAXX:1:3:1574:973/1
#bbaabbaaaabbaaaaaaaWaaaaaaZaaaaaaaaaaaaa_[_a_X_^UEEEE
                                          
if __name__ == '__main__':
    main()
