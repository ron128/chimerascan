'''
Created on May 13, 2011

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
import logging
import os

from chimerascan import pysam
from chimerascan.lib import config
from chimerascan.lib.sam import parse_pe_reads, parse_unpaired_pe_reads
from chimerascan.lib.seq import parse_fastq, fastq_to_string
from chimerascan.lib.transcriptome import build_gene_interval_trees, get_overlapping_genes
#from chimerascan.lib.transcriptome import build_exon_interval_trees, get_transcript_coords

def imin2(a, b):
    return a if a <= b else b

def get_min_genomic_distance(pe_reads):
    mindist = None
    for r1 in pe_reads[0]:
        for r2 in pe_reads[1]:
            if r1.rname != r2.rname:
                continue            
            d = imin2(abs(r1.aend - r2.pos), abs(r1.pos - r2.aend))
            if mindist is None:
                mindist = d
            else:
                mindist = imin2(mindist, d)
    return mindist

def get_concordant_hits(hitalns):
    # read may align to multiple locations ("hits")
    hit_tx_names = set(aln.g.tx_name for aln in hitalns[0])
    hit_gene_names = set(aln.g.gene_name for aln in hitalns[0])
    for segalns in hitalns[1:]:
        # each read hit may align with gaps to multiple segments
        if len(segalns) == 0:
            hit_tx_names = set()
            hit_gene_names = set()
        else:
            hit_tx_names.intersection_update(set(aln.g.tx_name for aln in segalns))
            hit_gene_names.intersection_update(set(aln.g.gene_name for aln in segalns))
        if len(hit_gene_names) == 0:
            break
    return hit_tx_names, hit_gene_names

def get_concordant_reads(readalns):
    if len(readalns) == 0:
        return set(),set()
    # process read1/read2
    read_tx_names, read_gene_names = set(), set()
    for hitalns in readalns:
        if len(hitalns) == 0:
            continue
        tx_names, gene_names = get_concordant_hits(hitalns)
        read_tx_names.update(tx_names)
        read_gene_names.update(gene_names)
    return read_tx_names, read_gene_names

def get_concordant_pairs(pe_alns):
    """
    return (tx_names, gene_names) tuple if there is a pair of reads that 
    aligns to the same transcript/gene in the correct orientation
    """
    r1_tx_names, r1_gene_names = get_concordant_reads(pe_alns[0])
    r2_tx_names, r2_gene_names = get_concordant_reads(pe_alns[1])
    return (r1_tx_names.intersection(r2_tx_names), 
            r1_gene_names.intersection(r2_gene_names))

def has_transcript_alignments(readalns):
    for hitalns in readalns:
        for segalns in hitalns:
            for aln in segalns:
                return True
    return False

def has_paired_transcript_alignments(pe_alns):
    return (has_transcript_alignments(pe_alns[0]) and
            has_transcript_alignments(pe_alns[1]))

def is_concordant(bamfh, pe_reads, exon_intervals, exon_trees,
                  max_fragment_length):
    # check genomic distance
    mindist = get_min_genomic_distance(pe_reads)
    if mindist <= max_fragment_length:
        return True, True
    # lookup alignments in transcript space
    pe_alns = [[],[]]
    for readnum,reads in enumerate(pe_reads):
        readhits = [get_transcript_coords(bamfh, read, exon_intervals, exon_trees) 
                    for read in reads]
        pe_alns[readnum] = readhits
    if not has_paired_transcript_alignments(pe_alns):
        # TODO: handle discordant unannotated candidates
        if mindist <= max_fragment_length:
            return True, True
        return True, True
    # screen for concordant reads
    tx_names, gene_names = get_concordant_pairs(pe_alns)
    # debugging output
#    r1_tx_names, r1_gene_names = get_concordant_reads(pe_alns[0])
#    r2_tx_names, r2_gene_names = get_concordant_reads(pe_alns[1])
#    if r1_gene_names.isdisjoint(r2_gene_names):
#        print 'READ1', r1_tx_names, r1_gene_names
#        print 'READ2', r2_tx_names, r2_gene_names
#        for readnum,readalns in enumerate(pe_alns):
#            print 'READ', readnum, "ALNS", readalns
#            for hitnum,hitalns in enumerate(readalns):
#                print 'HIT', hitnum, "ALNS", hitalns
#                for segnum,segalns in enumerate(hitalns):
#                    print 'SEG', segnum, "ALNS", segalns
#                    for i,aln in enumerate(segalns):
#                        print 'READ', readnum, 'HIT', hitnum, 'SEG', segnum, 'ALN', i, aln
#        for r1 in pe_reads[0]:
#            print 'READ1', bamfh.getrname(r1.rname), r1
#        for r2 in pe_reads[1]:
#            print 'READ2', bamfh.getrname(r2.rname), r2  
    return len(tx_names) > 0, len(gene_names) > 0


def get_bam_reads_qname(pe_reads):
    if (len(pe_reads[0]) == 0) and (len(pe_reads[1]) == 0):
        return None
    elif len(pe_reads[0]) == 0:        
        return pe_reads[1][0].qname
    return pe_reads[0][0].qname

def write_pe_fastq(fqreads, outfq, suffix):
    print >>outfq[0], fastq_to_string(fqreads[0], suffix="%s1" % (suffix))
    print >>outfq[1], fastq_to_string(fqreads[1], suffix="%s2" % (suffix))

def synchronize_bam_fastq(bam_pe_reads, fastq_iters, outfq, suffix):
    """
    parses through sorted FASTQ files until encountering a read with
    the same qname as current BAM alignment record.  writes all unaligned
    reads to the 'outfq' file descriptor.
    
    returns 'False' if the last read encountered was unaligned, 'True' 
    otherwise 
    """
    fqreads = [it.next() for it in fastq_iters]
    fq_qname = fqreads[0].qname
    bam_qname = get_bam_reads_qname(bam_pe_reads)
    # fastq and bam records must be the same
    while bam_qname != fq_qname:
        # this read must be unmapped, so write it to 
        # the new fastq file
        write_pe_fastq(fqreads, outfq, suffix)
        # get the next fastq record
        fqreads = [it.next() for it in fastq_iters]
        fq_qname = fqreads[0].qname
    # check for unmapped reads ("one-mappers")
    has_unmapped_read = False
    if len(bam_pe_reads[0]) == 0 or len(bam_pe_reads[1]) == 0:
        has_unmapped_read = True
    elif (any(r.is_unmapped for r in bam_pe_reads[0]) or
          any(r.is_unmapped for r in bam_pe_reads[1])):
        has_unmapped_read = True
    if has_unmapped_read:
        write_pe_fastq(fqreads, outfq, suffix)
        return True
    return False

def process_tophat_alignments(fastq_files, bam_file, gene_file,
                              max_fragment_length,
                              output_fastq_files, 
                              output_bam_file,
                              unpaired=False,
                              suffix="/"):
    # index genes 
    exon_intervals, exon_trees = build_exon_interval_trees(gene_file)
    # open input files
    bamfh = pysam.Samfile(bam_file, "rb")
    if unpaired:
        bam_iter = parse_unpaired_pe_reads(bamfh)
    else:
        bam_iter = parse_pe_reads(bamfh)
    fastq_iters = [parse_fastq(open(fq)) for fq in fastq_files]
    # open output files
    outfq = [open(fq, "w") for fq in output_fastq_files]
    outbamfh = pysam.Samfile(output_bam_file, "wb", template=bamfh)
    # iterate through fastq files and bam file
    try:
        while True:
            bam_pe_reads = bam_iter.next()
            # synchronize fastq and bam and write unmapped reads to a file
            is_unaligned = synchronize_bam_fastq(bam_pe_reads, fastq_iters, 
                                                 outfq, suffix)
            if is_unaligned:
                continue
            # if loop reaches this point then we have a paired-end
            # read where both pairs align.  now need to check if
            # the alignment is discordant
            tx_concordant, gene_concordant = \
                is_concordant(bamfh, bam_pe_reads, exon_intervals, 
                              exon_trees, max_fragment_length)
            if not gene_concordant:
                for r in bam_pe_reads[0]:
                    outbamfh.write(r)
                for r in bam_pe_reads[1]:
                    outbamfh.write(r)
    except StopIteration:
        pass
    # finish remaining fastq lines
    try:
        while True:
            fqreads = [it.next() for it in fastq_iters]
            print >>outfq[0], fastq_to_string(fqreads[0])
            print >>outfq[1], fastq_to_string(fqreads[1])
    except StopIteration:
        pass
    return config.JOB_SUCCESS

def main():
    from optparse import OptionParser
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = OptionParser("usage: %prog [options] <index> <read1.fq> "
                          "<read2.fq> <accepted_hits.bam> "
                          "<unaligned_1.fq> <unaligned_2.fq> "
                          "<discordant.bam>")
    parser.add_option('--max-fragment-length', dest="max_fragment_length", 
                      type="int", default=1000)
    parser.add_option('--unpaired', action="store_true", dest="unpaired",
                      default=False)
    parser.add_option('--suffix', dest="suffix", default="/")
    options, args = parser.parse_args()
    index_dir = args[0]    
    fastq_files = args[1:3]
    bam_file = args[3]
    unmapped_fastq_files = args[4:6]
    discordant_bam_file = args[6]
    gene_file = os.path.join(index_dir, config.GENE_FEATURE_FILE)
    process_tophat_alignments(fastq_files, bam_file, gene_file,
                              options.max_fragment_length,
                              unmapped_fastq_files,
                              discordant_bam_file,
                              unpaired=options.unpaired,
                              suffix=options.suffix)

if __name__ == '__main__': 
    main()
