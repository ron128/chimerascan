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
from chimerascan.lib.sam import CIGAR_S, parse_pe_reads, parse_unpaired_pe_reads
from chimerascan.lib.seq import DNA_reverse_complement, parse_fastq, fastq_to_string
from chimerascan.lib.transcriptome import build_gene_interval_trees, get_overlapping_genes

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

def get_read_genes(bamfh, reads, gene_trees):
    read_tx_names = set()
    read_gene_names = set()
    for read in reads:
        hits = get_overlapping_genes(bamfh, read, gene_trees)
        hit_tx_names = None
        hit_gene_names = None
        for interval_genes in hits:
            if hit_gene_names is None:
                hit_tx_names = set(g.tx_name for g in interval_genes)
                hit_gene_names = set(g.gene_name for g in interval_genes)
            else:
                hit_tx_names.intersection_update(g.tx_name for g in interval_genes)
                hit_gene_names.intersection_update(g.gene_name for g in interval_genes)
        if hit_gene_names is not None:
            read_tx_names.update(hit_tx_names)
            read_gene_names.update(hit_gene_names)
    return read_tx_names, read_gene_names

def check_concordant(bamfh, pe_reads, gene_trees, max_fragment_length):
    # check genomic distance
    # TODO: account for introns in genomic distance checks
    # get overlapping genes
    mindist = get_min_genomic_distance(pe_reads)
    if mindist <= max_fragment_length:
        return True, True    
    r1_tx_names, r1_gene_names = \
        get_read_genes(bamfh, pe_reads[0], gene_trees)
    r2_tx_names, r2_gene_names = \
        get_read_genes(bamfh, pe_reads[1], gene_trees)
    pe_tx_names = r1_tx_names.intersection(r2_tx_names) 
    pe_gene_names = r1_gene_names.intersection(r2_gene_names)    
    return len(pe_tx_names) > 0, len(pe_gene_names) > 0

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
        return True, fqreads
    return False, fqreads

def write_discordant_reads(bamfh, fqreads, bam_pe_reads):
    for readnum in xrange(len(fqreads)):
        fq = fqreads[readnum]
        reads = bam_pe_reads[readnum]
        num_alignments = len(reads)
        for i,r in enumerate(reads):        
            # make sequence soft clipped
            ext_length = len(fq.seq) - len(r.seq)
            cigar_softclip = [(CIGAR_S, ext_length)]
            cigar = r.cigar
            # reconstitute full length sequence in read
            if r.is_reverse:
                #print 'REV before', r.qname, r.seq, r.qual, bamfh.getrname(r.rname), r.pos
                seq = DNA_reverse_complement(fq.seq)
                qual = fq.qual[::-1]
                r.seq = seq
                r.qual = qual
                if ext_length > 0:
                    r.cigar = cigar_softclip + cigar
                else:
                    r.cigar = cigar
                #print 'REV after', r.qname, r.seq, r.qual, bamfh.getrname(r.rname), r.pos
            else:
                #print 'FWD before', r.qname, r.seq, r.qual, bamfh.getrname(r.rname), r.pos
                r.seq = fq.seq
                r.qual = fq.qual
                if ext_length > 0:
                    r.cigar = cigar + cigar_softclip
                else:
                    r.cigar = cigar
                #print 'FWD after', r.qname, r.seq, r.qual, bamfh.getrname(r.rname), r.pos
            tagdict = dict(r.tags)
            # TODO: bug in pysam handling CP tag, fix by forcing to integer
            if "CP" in tagdict:
                tagdict["CP"] = int(tagdict["CP"])
            # annotate reads with 'HI', and 'IH' tags
            r.tags = [("HI",i), ("IH",num_alignments)] + tagdict.items()
            bamfh.write(r)


def find_discordant_tophat_alignments(fastq_files, bam_file, gene_file,
                                      max_fragment_length,
                                      output_fastq_files, 
                                      output_bam_file,
                                      unpaired=False,
                                      suffix="/"):
    # index genes 
    gene_trees = build_gene_interval_trees(gene_file)
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
            is_unaligned,fqreads = synchronize_bam_fastq(bam_pe_reads, fastq_iters, 
                                                         outfq, suffix)
            if is_unaligned:
                continue
            # if loop reaches this point then we have a paired-end
            # read where both pairs align.  now need to check if
            # the alignment is discordant
            tx_concordant, gene_concordant = \
                check_concordant(bamfh, bam_pe_reads, gene_trees, 
                                 max_fragment_length)
            if not gene_concordant:
                write_discordant_reads(outbamfh, fqreads, bam_pe_reads)
    except StopIteration:
        pass
    # finish remaining fastq lines
    try:
        while True:
            fqreads = [it.next() for it in fastq_iters]
            write_pe_fastq(fqreads, outfq, suffix)
    except StopIteration:
        pass
    bamfh.close()
    for fh in outfq:
        fh.close()
    outbamfh.close()
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
    find_discordant_tophat_alignments(fastq_files, bam_file, gene_file,
                                      options.max_fragment_length,
                                      unmapped_fastq_files,
                                      discordant_bam_file,
                                      unpaired=options.unpaired,
                                      suffix=options.suffix)

if __name__ == '__main__': 
    main()
