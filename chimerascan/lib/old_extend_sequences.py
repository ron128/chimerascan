'''
Created on Jan 23, 2011

@author: mkiyer
'''
import logging
import re

# local imports
from find_discordant_reads import Chimera

def parse_qname_file(qname_fh):
    for line in qname_fh:
        qname, mate = line.strip().split('\t')
        yield qname, int(mate)

def parse_discordant_by_qname(infh):
    chimeras = [] 
    for line in infh:
        chimera = Chimera.from_bedpe(line)
        qname = chimera.qname
        if len(chimeras) > 0 and (qname != chimeras[-1].qname):
            yield chimeras
            chimeras = []
        chimeras.append(chimera)
    if len(chimeras) > 0:
        yield chimeras

def parse_fastq(line_iter):
    try:        
        qname = line_iter.next().rstrip()[1:]
        newqname = re.split(r'/\d$', qname)[0]
        suffix_length = len(qname) - len(newqname)                    
        seq = line_iter.next().rstrip()
        line_iter.next()
        qual = line_iter.next().rstrip()
        yield newqname, seq, qual
        while True:
            # qname
            qname = line_iter.next().rstrip()[1:]
            qname = qname[:len(qname)-suffix_length]
            # seq
            seq = line_iter.next().rstrip()
            # qname again (skip)
            line_iter.next()
            # qual
            qual = line_iter.next().rstrip()
            yield qname, seq, qual
    except StopIteration:
        pass
    
def fetch_fastq_qname(fastq_iters, qname):
    while True:
        records = [it.next() for it in fastq_iters]
        if records[0][0] == qname:
            break
    return records

def to_fastq(mate, qname, seq, qual):
    return "@%s/%d\n%s\n+%s/%d\n%s" % (qname, mate, seq, qname, mate, qual)

def extend_sequences(input_fastq_files, discordant_bedpe_file,
                     output_discordant_bedpe_file, 
                     output_spanning_fastq_file):
    fastq_iters = [parse_fastq(open(fq)) for fq in input_fastq_files]
    bedpe_fh = open(output_discordant_bedpe_file, "w")
    fastq_fh = open(output_spanning_fastq_file, "w")        
    num_lines = 0
    for chimeras in parse_discordant_by_qname(open(discordant_bedpe_file)):        
        # find sequences in fastq file
        fastq_records = fetch_fastq_qname(fastq_iters, chimeras[0].qname)
        print >>fastq_fh, to_fastq(1, *fastq_records[0]) 
        print >>fastq_fh, to_fastq(2, *fastq_records[1])
        for chimera in chimeras:
            mate1, mate2 = (0, 1) if chimera.read1_is_5prime else (1, 0)
            # get extended sequences
            seq1, seq2 = fastq_records[mate1][1], fastq_records[mate2][1]
            qual1, qual2 = fastq_records[mate1][2], fastq_records[mate2][2]
            # add extended sequences to chimera
            chimera.mate5p.seq = seq1
            chimera.mate5p.qual = qual1
            chimera.mate3p.seq = seq2
            chimera.mate3p.qual = qual2  
            print >>bedpe_fh, chimera.to_bedpe()
        num_lines += 1
    fastq_fh.close()
    bedpe_fh.close()
    logging.info("[EXTENDSEQ] Processed %d discordant reads" % (num_lines))
    

def main():
    from optparse import OptionParser
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = OptionParser("usage: %prog [options]")
    options, args = parser.parse_args()    
    discordant_bedpe_file = args[0]
    output_discordant_bedpe_file = args[1]
    output_spanning_fastq_file = args[2]
    input_fastq_files = args[3:]
    extend_sequences(input_fastq_files, discordant_bedpe_file,
                     output_discordant_bedpe_file, 
                     output_spanning_fastq_file)

if __name__ == '__main__':
    main()

