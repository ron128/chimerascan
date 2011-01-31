'''
Created on Jan 29, 2011

@author: mkiyer
'''
import logging
import re

# local imports
from find_discordant_reads import DiscordantFragment

def parse_discordant_by_qname(infh):
    frags = [] 
    for line in infh:
        fields = line.strip().split('\t')
        frag = DiscordantFragment.from_list(fields)        
        qname = frag.qname
        if len(frags) > 0 and (qname != frags[-1].qname):
            yield frags
            frags = []
        frags.append(frag)
    if len(frags) > 0:
        yield frags

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

def extend_sequences(input_fastq_files, 
                     discordant_bedpe_file,
                     output_discordant_bedpe_file): 
    fastq_iters = [parse_fastq(open(fq)) for fq in input_fastq_files]
    bedpe_fh = open(output_discordant_bedpe_file, "w")
    num_lines = 0
    for discord_reads in parse_discordant_by_qname(open(discordant_bedpe_file)):        
        # find sequences in fastq file
        fastq_records = fetch_fastq_qname(fastq_iters, discord_reads[0].qname)
        for discord_read in discord_reads:
            mate5p, mate3p = (0, 1) if discord_read.read1_is_sense else (1, 0)
            # get extended sequences
            seq5p, seq3p = fastq_records[mate5p][1], fastq_records[mate3p][1]
            qual5p, qual3p = fastq_records[mate5p][2], fastq_records[mate3p][2]            
            # add extended sequences to chimera
            discord_read.clust5p.seq = seq5p
            discord_read.clust5p.qual = qual5p
            discord_read.clust3p.seq = seq3p
            discord_read.clust3p.qual = qual3p
            print >>bedpe_fh, '\t'.join(map(str, discord_read.to_list()))
        num_lines += 1
    bedpe_fh.close()
    logging.info("[EXTENDSEQ] Processed %d discordant reads" % (num_lines))
    

def main():
    from optparse import OptionParser
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = OptionParser("usage: %prog [options] <input.bedpe> <output.bedpe> <fastq_file(s)")
    options, args = parser.parse_args()    
    discordant_bedpe_file = args[0]
    output_discordant_bedpe_file = args[1]
    input_fastq_files = args[2:]
    extend_sequences(input_fastq_files, 
                     discordant_bedpe_file,
                     output_discordant_bedpe_file)

if __name__ == '__main__':
    main()

