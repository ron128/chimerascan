'''
Created on Jan 30, 2011

@author: mkiyer
'''
import collections
import logging

from nominate_chimeras2 import Chimera, parse_discordant_reads
from find_discordant_reads2 import DiscordantType

def to_fastq(mate, qname, seq, qual):
    return "@%s/%d\n%s\n+%s/%d\n%s" % (qname, mate+1, seq, qname, mate+1, qual)

def is_spanning(start, end, juncs):
    return any(start < j < end for j in juncs)

def check_fragment(frag, tx5p, tx3p):
    # mates that have two split segments mapping discordantly
    # are automatically considered for spanning read detection
    write5p = frag.discordant_type.discordant5p
    write3p = frag.discordant_type.discordant3p
    # check the padding boundaries to find mates with positions 
    # past the junction
    if (not write5p) and (frag.clust5p.rname in tx5p):
        write5p = is_spanning(frag.clust5p.pad_start, 
                              frag.clust5p.pad_end,
                              tx5p[frag.clust5p.rname])
    if (not write3p) and (frag.clust3p.rname in tx3p):
        write3p = is_spanning(frag.clust3p.pad_start, 
                              frag.clust3p.pad_end,
                              tx5p[frag.clust3p.rname])
    # check "one-mapper" fragments
    if frag.discordant_type.code == DiscordantType.CONCORDANT_SINGLE:
        # one of mates mapped and other is unmapped, so check
        # the mapped mate and see whether it matches a chimera
        # candidate 
        # TODO: check junction position to further refine decision
        # by omitting reads that are far from the predicted junction
        if (frag.clust5p.rname == "*") and (frag.clust3p.rname in tx3p):
            print 'CONCORDANT_SINGLE', frag.qname
            write5p = True
        if (frag.clust3p.rname == "*") and (frag.clust5p.rname in tx5p):
            print 'CONCORDANT_SINGLE', frag.qname
            write3p = True
    # write the potential spanning reads
    mate = int(frag.read1_is_sense)
    read1, read2 = None, None
    if write5p:
        fq = to_fastq(mate, frag.qname, frag.clust5p.seq, frag.clust5p.qual)            
        if frag.read1_is_sense:
            read1 = fq 
        else:
            read2 = fq
    if write3p:
        fq = to_fastq(mate, frag.qname, frag.clust3p.seq, frag.clust3p.qual)            
        if frag.read1_is_sense:
            read2 = fq 
        else:
            read1 = fq
    return read1, read2

def nominate_spanning_reads(discordant_reads_fh,
                            chimeras_fh,
                            fastq_fh):
    # build index of chimera candidates
    logging.info("Indexing chimera candidates")
    tx5p = collections.defaultdict(lambda: [])
    tx3p = collections.defaultdict(lambda: [])
    for chimera in Chimera.parse(chimeras_fh):
        tx5p[chimera.mate5p.tx_name].append(chimera.mate5p.junc_pos)
        tx3p[chimera.mate3p.tx_name].append(chimera.mate3p.junc_pos)
    # parse discordant reads    
    logging.info("Nominating spanning reads")    
    read1, read2 = None, None
    prev_qname = None
    for frag in parse_discordant_reads(discordant_reads_fh):        
        if frag.discordant_type.is_genome:
            continue
        qname = frag.qname
        if prev_qname is not None and (qname != prev_qname):
            if read1 is not None:
                print >>fastq_fh, read1
            if read2 is not None:
                print >>fastq_fh, read2
            read1, read2 = None, None
        prev_qname = qname
        if (read1 is not None) and (read2 is not None):
            continue
        r1, r2 = check_fragment(frag, tx5p, tx3p)
        if read1 is None: read1 = r1
        if read2 is None: read2 = r2
    if read1 is not None:
        print >>fastq_fh, read1
    if read2 is not None:
        print >>fastq_fh, read2


def main():
    from optparse import OptionParser
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = OptionParser("usage: %prog [options] <qname_sorted_discordant_reads> <chimeras> <output.fq>")
    options, args = parser.parse_args()
    discordant_reads_file = args[0]
    chimeras_file = args[1]
    fastq_file = args[2]
    nominate_spanning_reads(open(discordant_reads_file, 'r'),
                            open(chimeras_file, 'r'),
                            open(fastq_file, 'w'))

if __name__ == '__main__':
    main()
