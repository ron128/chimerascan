'''
Created on Jan 30, 2011

@author: mkiyer

chimerascan: chimeric transcript discovery using RNA-seq

Copyright (C) 2011 Matthew Iyer

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''
import collections
import logging
from bisect import bisect

from chimerascan import pysam
from chimerascan.lib import config
from chimerascan.lib.sam import parse_pe_reads
from chimerascan.lib.chimera import Chimera, DiscordantRead

def to_fastq(qname, readnum, seq, qual):
    return "@%s/%d\n%s\n+\n%s" % (qname, readnum+1, seq, qual)

def is_spanning(start, end, juncs):
    return any(start < j < end for j in juncs)

def check_fragment(frag, tx5p, tx3p, nonmapping=True):
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
    if nonmapping and (frag.discordant_type.code == DiscordantType.NONMAPPING):
        # TODO: automatically completely non-mapping reads that may
        # be double-overlapping spanning reads, but only do this in
        # single-segment mode to increase sensitivity
        write5p = True
        write3p = True
    elif frag.discordant_type.code == DiscordantType.CONCORDANT_SINGLE:
        # one of mates mapped and other is unmapped, so check
        # the mapped mate and see whether it matches a chimera
        # candidate 
        # TODO: check junction position to further refine decision
        # by omitting reads that are far from the predicted junction
        if (frag.clust5p.rname == "*") and (frag.clust3p.rname in tx3p):
            write5p = True
        if (frag.clust3p.rname == "*") and (frag.clust5p.rname in tx5p):
            write3p = True
    # write the potential spanning reads
    reads = [None, None]
    if write5p:
        mate = 0 if frag.read1_is_sense else 1
        reads[mate] = to_fastq(mate, frag.qname, frag.clust5p.seq, frag.clust5p.qual)
    if write3p:
        mate = 1 if frag.read1_is_sense else 0
        reads[mate] = to_fastq(mate, frag.qname, frag.clust3p.seq, frag.clust3p.qual)
    return reads

def nominate_spanning_reads2(discordant_reads_fh,
                            chimeras_fh,
                            fastq_fh):
    # build index of chimera candidates
    logging.info("Indexing chimera candidates")
    tx5p = collections.defaultdict(lambda: [])
    tx3p = collections.defaultdict(lambda: [])
    for chimera in Chimera.parse(chimeras_fh):
        tx5p[chimera.mate5p.tx_name].append(chimera.mate5p.end)
        tx3p[chimera.mate3p.tx_name].append(chimera.mate3p.start)
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
        # skip if reads already found
        if (read1 is not None) and (read2 is not None):
            continue
        # update read fastq
        r1, r2 = check_fragment(frag, tx5p, tx3p)
        if read1 is None: read1 = r1
        if read2 is None: read2 = r2
        prev_qname = qname
    if read1 is not None:
        print >>fastq_fh, read1
    if read2 is not None:
        print >>fastq_fh, read2

def nominate_spanning_reads(chimera_file, unmapped_bam_file, output_fastq_file):
    # find all reads that need to be remapped to see if they span the 
    # breakpoint junction
    fqfh = open(output_fastq_file, "w")
    remap_qnames = set()
    breaks5p = collections.defaultdict(lambda: [])
    breaks3p = collections.defaultdict(lambda: [])
    for c in Chimera.parse(open(chimera_file)):
        end5p = c.partner5p.end
        start3p = c.partner3p.start
        # keep track of all breakpoints
        breaks5p[c.partner5p.tx_name].append(end5p)
        breaks3p[c.partner5p.tx_name].append(start3p)
        for r5p,r3p in c.encomp_read_pairs:            
            # if 5' read overlaps breakpoint then it should be remapped
            if r5p.clipstart < end5p < r5p.clipend:
                key5p = (r5p.qname, r5p.readnum)
                if key5p not in remap_qnames:
                    remap_qnames.add((r5p.qname, r5p.readnum))
                    print >>fqfh, to_fastq(r5p.qname, r5p.readnum, 
                                           r5p.seq, "I" * len(r5p.seq))
            # if 3' read overlaps breakpoint then it should be remapped
            if r3p.clipstart < start3p < r3p.clipend:
                key3p = (r3p.qname, r3p.readnum)
                if key3p not in remap_qnames:
                    remap_qnames.add((r3p.qname, r3p.readnum))
                    print >>fqfh, to_fastq(r3p.qname, r3p.readnum, 
                                           r3p.seq, "I" * len(r3p.seq))
    # sort breakpoint positions within each gene
    for tx_name in breaks5p.keys():
        breaks5p[tx_name] = sorted(breaks5p[tx_name])
    for tx_name in breaks3p.keys():
        breaks3p[tx_name] = sorted(breaks3p[tx_name])   
    # check read pairs with one or both unmapped, and remap those 
    # as well
    bamfh = pysam.Samfile(unmapped_bam_file, "rb")
    for pe_reads in parse_pe_reads(bamfh):
        for readnum in xrange(0, 2):
            print >>fqfh, to_fastq(pe_reads[readnum][0].qname, readnum, 
                                   pe_reads[readnum][0].seq,
                                   pe_reads[readnum][0].qual) 
#            # add unmapped reads
#            if reads[0].is_unmapped:
#                readnum = 2 if reads[0].is_read2 else 1
#                print >>fqfh, to_fastq(reads[0].qname, readnum, reads[0].seq, 
#                                       "I" * len(reads[0].seq))
#                # TODO: remove this
#                assert len(reads) == 1
#            else:                
#                remap = False
#                for r in reads:
#                    tx_name = config.GENE_REF_PREFIX + bamfh.getrname(r.rname)
#                    # check if this read overlaps a breakpoint
#                    
#                    bisect()
    bamfh.close()
    return config.JOB_SUCCESS

def main():
    from optparse import OptionParser
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = OptionParser("usage: %prog [options] <chimeras.txt> <unmapped_reads.bam> <out.fq>")
    options, args = parser.parse_args()
    chimera_file = args[0]
    bam_file = args[1]
    output_fastq_file = args[2]
    nominate_spanning_reads(chimera_file, bam_file, output_fastq_file)

if __name__ == '__main__':
    main()
