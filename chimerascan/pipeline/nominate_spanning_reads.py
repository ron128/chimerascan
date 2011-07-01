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

from chimerascan import pysam
from chimerascan.lib import config
from chimerascan.lib.sam import parse_pe_reads
from chimerascan.lib.chimera import Chimera, DiscordantRead
from chimerascan.lib.config import GENE_REF_PREFIX

def to_fastq(qname, readnum, seq, qual):
    return "@%s/%d\n%s\n+\n%s" % (qname, readnum+1, seq, qual)

def nominate_encomp_spanning_reads(chimera_file, output_fastq_file):
    # find all reads that need to be remapped to see if they span the 
    # breakpoint junction
    fqfh = open(output_fastq_file, "w")
    remap_qnames = set()
    #breaks5p = collections.defaultdict(lambda: [])
    #breaks3p = collections.defaultdict(lambda: [])
    for c in Chimera.parse(open(chimera_file)):
        end5p = c.partner5p.end
        start3p = c.partner3p.start
        # keep track of all breakpoints
        #breaks5p[GENE_REF_PREFIX + c.partner5p.tx_name].append(end5p)
        #breaks3p[GENE_REF_PREFIX + c.partner3p.tx_name].append(start3p)
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
    #for rname in breaks5p.keys():
    #    breaks5p[rname] = sorted(breaks5p[rname])
    #for tx_name in breaks3p.keys():
    #    breaks3p[rname] = sorted(breaks3p[rname])
    fqfh.close()
    return config.JOB_SUCCESS

def nominate_unmapped_spanning_reads(unmapped_bam_file, output_fastq_file): 
    # find all reads that need to be remapped to see if they span the 
    # breakpoint junction
    fqfh = open(output_fastq_file, "w")
    # check read pairs with one or both unmapped, and remap those 
    # as well
    bamfh = pysam.Samfile(unmapped_bam_file, "rb")    
    for pe_reads in parse_pe_reads(bamfh):
        # remap all unmapped reads
        for readnum,reads in enumerate(pe_reads):
            if any(r.is_unmapped for r in reads):
                print >>fqfh, to_fastq(pe_reads[readnum][0].qname, readnum, 
                                       pe_reads[readnum][0].seq,
                                       pe_reads[readnum][0].qual) 

    bamfh.close()
    fqfh.close()
    return config.JOB_SUCCESS

def main():
    from optparse import OptionParser
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = OptionParser("usage: %prog [options] <chimeras.txt> "
                          "<unmapped_reads.bam> <encomp_remap.fq> "
                          "<unmapped_remap.fq>")
    options, args = parser.parse_args()
    chimera_file = args[0]
    bam_file = args[1]
    encomp_remap_fastq_file = args[2]
    spanning_fastq_file = args[3]
    nominate_encomp_spanning_reads(chimera_file, encomp_remap_fastq_file)
    nominate_unmapped_spanning_reads(bam_file, spanning_fastq_file) 

if __name__ == '__main__':
    main()
