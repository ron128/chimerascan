'''
Created on Oct 22, 2010

@author: mkiyer

nohup time python /exds/users/mkiyer/sequel/sw/sequel2/tools/rnaseq/fusionpipeline/align.py --bowtie-bin /exds/users/mkiyer/sequel/sw/alignment/bowtie/bowtie-0.12.7/bowtie --bowtie-index /exds/users/mkiyer/sequel/alignment_indexes/bowtie/hg19_genome_knowngene_combined/hg19_genome_knowngene_combined --gene-bed /exds/users/mkiyer/sequel/reference_sequences/hg19/annotations/knownGene_2010_Oct_14.bed --gene-fasta-prefix hg19_knownGene_ --multihits 100 --mismatches 2 --seed-length 25 --quals solexa1.3-quals --bowtie-threads 3 s_3_1_sequence.txt s_3_2_sequence.txt discordant_reads.bam > align.log 2>&1 &

'''
import argparse
import tempfile
import os
import sys
import logging
import collections
import itertools
import subprocess
import multiprocessing
import shutil
import numpy as np

from bx.intervals.intersection import Interval, IntervalTree
import pysam

from base import parse_bed12_file
from gene_to_genome import build_gene_to_genome_map, transcriptome_to_genome

class ExonData(object):
    __slots__ = ('id', 'start', 'end', 'strand', 'reads', 'ambiguous_reads', 'weighted_cov', 'weighted_ambiguous_cov')
    def __init__(self):
        self.reads = 0
        self.ambiguous_reads = 0
        self.weighted_cov = 0.0
        self.weighted_ambiguous_cov = 0.0

def build_interval_trees(samfh, bedfile):
    # add exons to interval tree and data array    
    gene_trees = collections.defaultdict(lambda: IntervalTree())    
    current_exon_id = 0
    exon_trees = collections.defaultdict(lambda: IntervalTree())
    exon_keys = collections.defaultdict(lambda: {})
    exon_data = []
    exon_gene_map = collections.defaultdict(lambda: [])
    rname_tid_map = dict((rname,i) for i,rname in enumerate(samfh.references))    
    for g in parse_bed12_file(open(bedfile)):
        if g.chrom not in rname_tid_map:
            continue
        # get reference index in sam file
        tid = rname_tid_map[g.chrom]
        # add gene to interval tree
        gene_interval = Interval(g.tx_start, g.tx_end, strand=g.strand, value=g.name)
        gene_trees[tid].insert_interval(gene_interval)
        # add exons to interval tree
        for exon_num, exon in enumerate(g.exons):
            start, end = exon
            key = (start, end, g.strand)            
            if key not in exon_keys[tid]:
                exon_id = current_exon_id
                current_exon_id += 1                
                exon_keys[tid][key] = exon_id
                eobj = ExonData()
                eobj.id = exon_id
                eobj.start = start
                eobj.end = end
                eobj.strand = g.strand
                exon_trees[tid].insert_interval(Interval(start, end, strand=g.strand, value=eobj))
                exon_data.append(eobj)
            else:
                exon_id = exon_keys[tid][key]
            # store lookup between exons to genes
            exon_gene_map[exon_id].append((g.name, exon_num))
    return gene_trees, exon_trees, exon_data, exon_gene_map

def make_fifo(base_dir):
    tmpdir = tempfile.mkdtemp(suffix='fifo', prefix='tmp', dir=base_dir)
    fifo_file = os.path.join(tmpdir, "fifo")
    try:
        os.mkfifo(fifo_file)
    except OSError, e:
        logging.error("Failed to create FIFO: %s" % e)
    return tmpdir, fifo_file

def make_temp(base_dir, suffix=''):
    fd,name = tempfile.mkstemp(suffix=suffix, prefix='tmp', dir=base_dir)
    os.close(fd)
    return name

def iter_multihit_alignments(samfh):    
    reads = []
    for read in samfh:
        if len(reads) > 0 and read.qname != reads[-1].qname:
            yield reads
            reads = []
        reads.append(read)
    if len(reads) > 0:
        yield reads

CIGAR_M = 0 #match  Alignment match (can be a sequence match or mismatch)
CIGAR_I = 1 #insertion  Insertion to the reference
CIGAR_D = 2 #deletion  Deletion from the reference
CIGAR_N = 3 #skip  Skipped region from the reference
CIGAR_S = 4 #softclip  Soft clip on the read (clipped sequence present in <seq>)
CIGAR_H = 5 #hardclip  Hard clip on the read (clipped sequence NOT present in <seq>)
CIGAR_P = 6 #padding  Padding (silent deletion from the padded reference sequence)

def get_aligned_read_intervals(read):
    intervals = []
    # insert read into cluster tree
    astart,aend = read.pos, read.pos
    for op,length in read.cigar:
        if length == 0: continue
        if (op == CIGAR_I) or (op == CIGAR_S) or (op == CIGAR_H): continue
        if (op == CIGAR_P): assert False 
        if (op == CIGAR_N):
            assert astart != aend
            intervals.append((astart, aend))
            #print read.qname, read.cigar, ref, astart, aend
            astart = aend + length
        aend += length
    assert astart != aend
    if aend > astart:
        #print read.qname, read.cigar, ref, astart, aend
        intervals.append((astart, aend))
    assert aend == read.aend
    return intervals

#def get_genes_overlapping_reads(reads, gene_trees):
#    hits = set()
#    for read in reads:
#        if read.is_unmapped or read.is_qcfail:
#            continue
#        rname = read.rname
#        for interval in get_aligned_read_intervals(read):            
#            hits.update([hit.value for hit in gene_trees[rname].find(interval[0], interval[1])])
#    return hits

def assign_coverage(reads, gene_trees, exon_trees):
    # count mapped reads
    mapped_reads = [r for r in reads
                    if not (r.is_unmapped or r.is_qcfail)]
    # use equal weighting to assign coverage for multimapping reads
    weighted_cov = 0
    if len(mapped_reads) > 0:
        weighted_cov = 1.0 / float(len(mapped_reads))
    # intersect each mapped read with exon interval trees
    # and also find gene names
    gene_names = set()
    for read in mapped_reads:
        rname = read.rname
        # get the distinct genomic intervals that the
        # read aligns to
        intervals = get_aligned_read_intervals(read)
        interval_exon_sets = []
        for start,end in intervals:
            # update gene names
            gene_names.update(hit.value for hit in gene_trees[rname].find(start, end))
            # find exons that contain this read
            eobjs = [hit.value for hit in exon_trees[rname].find(start, end)
                     if start >= hit.start and end <= hit.end]
            interval_exon_sets.append(eobjs)        
        for interval_exons in interval_exon_sets:
            if len(interval_exons) == 0:
                continue
            ambiguous = len(interval_exons) > 1
            if ambiguous:
                cov = weighted_cov / (len(interval_exons) * len(interval_exon_sets))
            else:
                cov = weighted_cov / (len(interval_exons))
            for eobj in interval_exons:
                if ambiguous:        
                    eobj.ambiguous_reads += 1
                    eobj.weighted_ambiguous_cov += cov
                else:
                    eobj.reads += 1
                    eobj.weighted_cov += cov 
    return gene_names, len(mapped_reads)

def check_reads_within_insert_size_range(mate1_reads, mate2_reads, insert_size_max):
    # bin reads by reference name
    read1_by_rname = collections.defaultdict(lambda: set())
    read2_by_rname = collections.defaultdict(lambda: set())
    for read1 in mate1_reads:
        if read1.is_unmapped or read1.is_qcfail:
            return False
        read1_by_rname[read1.rname].add(read1)
    for read2 in mate2_reads:
        if read2.is_unmapped or read2.is_qcfail:
            return False
        read2_by_rname[read2.rname].add(read2)
    # find the smallest insert size between any two reads
    min_isize = insert_size_max + 1    
    # search by reference and find fusion reads
    for rname,read1_set in read1_by_rname.iteritems():
        # see if there is a match on the same chromosome
        if rname in read2_by_rname:
            # check for matches to the same gene
            read2_set = read2_by_rname[rname]
            for r1 in read1_set:
                for r2 in read2_set:
                    if r1.pos > r2.pos:
                        isize = r1.pos - r2.aend
                    else:
                        isize = r2.pos - r1.aend
                    min_isize = min(min_isize, isize)
                    if min_isize <= insert_size_max:
                        return True            
    return False

def set_read_flags(r, mate, num_hits):
    if mate == 0:
        r.is_read1 = True
    else:
        r.is_read2 = True
    r.tags = r.tags + [('NH', num_hits)]
    r.qname = r.qname[:-2]
    r.is_paired = True
    r.is_proper_pair = False
    r.mate_is_unmapped = True
    r.mrnm = -1
    r.mpos = -1
    r.isize = 0

def parse_sam_files(samfhs, maxlen=100000):
    buf = []
    for i in xrange(maxlen):
        buf.append(([],[]))
    buf_size = 0
    next_buf_ind = 0
    buf_ind = 0
    qname_ind_map = {}
    sam_iters = [iter_multihit_alignments(fh) for fh in samfhs]
    #for read_segments in itertools.izip(*sam_iters):
    try:
        while True:
            read_pair = [sam_iter.next() for sam_iter in sam_iters]
            for mate,reads in enumerate(read_pair):
                for read in reads:
                    # set read flags
                    if mate == 0:
                        read.is_read1 = True
                    else:
                        read.is_read2 = True
                    read.qname = read.qname[:-2]
                    read.is_paired = True
                    read.is_proper_pair = False
                    read.mate_is_unmapped = True
                    read.mrnm = -1
                    read.mpos = -1
                    read.isize = 0
                    # add read to buffer
                    if read.qname not in qname_ind_map:
                        if buf_size == maxlen:
                            return_read_pair = buf[next_buf_ind]
                            del qname_ind_map[return_read_pair[0][0].qname]
                            yield buf[next_buf_ind]
                        else:
                            buf_size += 1
                        buf_ind = next_buf_ind
                        next_buf_ind += 1
                        if next_buf_ind == maxlen:
                            next_buf_ind = 0
                        qname_ind_map[read.qname] = buf_ind
                        buf[buf_ind] = ([],[])
                    else:
                        buf_ind = qname_ind_map[read.qname]
                    buf[buf_ind][mate].append(read)
    except StopIteration:
        pass
    for buf_ind in xrange(buf_size):
        yield buf[buf_ind]

def write_discordant_reads(reads, num_hits, outfh):
    multihit = False
    for r in reads:
        # do not want more than one unmapped read per mate
        if r.is_unmapped and multihit:
            continue
        multihit = True
        #logging.debug("read: %s num_hits: %d" % (r.qname, num_hits))
        r.tags = r.tags + [('NH', num_hits)]
        outfh.write(r)

def write_exon_expr_data(exon_data, exon_gene_map, outfh):
    # need to write the exon coverage data to a file
    print >>outfh, '\t'.join(["#gene_name", "exon_number", 
                              "exon_start", "exon_end", "exon_strand",
                              "num_reads", "num_ambiguous_reads",
                              "weighted_cov", "weighted_ambiguous_cov"])
    for eobj in exon_data:
        for gene_name, exon_num in exon_gene_map[eobj.id]:
            print >>outfh, '\t'.join(map(str, [gene_name,
                                               exon_num,
                                               eobj.start,
                                               eobj.end,
                                               eobj.strand,
                                               eobj.reads,
                                               eobj.ambiguous_reads,
                                               eobj.weighted_cov,
                                               eobj.weighted_ambiguous_cov]))

def join_paired_end_reads(sam_files, output_sam_file, output_expr_file,
                          gene_bed_file, insert_size_max):
    debug_count = 0
    debug_every = 1e5
    debug_next = debug_every
    try:
        samfhs = [pysam.Samfile(f, "r") for f in sam_files]
        header = samfhs[0].header
        outfh = pysam.Samfile(output_sam_file, "wb", template=samfhs[0])
        # build interval trees for comparing hits to the same gene
        logging.info("Constructing interval intersection indexes...")
        gene_trees, exon_trees, exon_data, exon_gene_map = build_interval_trees(outfh, gene_bed_file)
        # iterate through paired-end alignments
        logging.info("Processing paired alignments...")
        for mate1_reads, mate2_reads in parse_sam_files(samfhs):
            debug_count += 1
            if debug_count == debug_next:
                debug_next += debug_every
                logging.info("Processed %d reads" % debug_count)
#            for r in mate1_reads:
#                logging.info("mate1 %s" % (r.qname))
#            for r in mate2_reads:
#                logging.info("mate2 %s" % (r.qname))
            # find all gene hits for each mate
#            r1_genes = get_genes_overlapping_reads(mate1_reads, gene_trees)
#            r2_genes = get_genes_overlapping_reads(mate2_reads, gene_trees)
            # assign coverage
            r1_genes, mate1_numhits = assign_coverage(mate1_reads, gene_trees, exon_trees)
            r2_genes, mate2_numhits = assign_coverage(mate2_reads, gene_trees, exon_trees)
            # intersect gene lists
            shared_genes = r1_genes.intersection(r2_genes)
            # if there are shared genes, then we exclude this mate pair
            # from being a gene fusion candidate
            if len(shared_genes) > 0:
                continue
            # if the insert size between the mates is small, then exclude
            # from being a gene fusion candidate
            isize_filter = check_reads_within_insert_size_range(mate1_reads, mate2_reads, 
                                                                insert_size_max)
            if isize_filter:
                continue
            # write fusion reads with all the mate1 reads followed
            # by all the mate2 reads.  the reads are still not
            # paired at this point since this could result in a 
            # combinatorial explosion
            write_discordant_reads(mate1_reads, mate1_numhits, outfh)
            write_discordant_reads(mate2_reads, mate2_numhits, outfh)
    except StopIteration:
        pass
    write_exon_expr_data(exon_data, exon_gene_map, open(output_expr_file, "w"))
    

def setup_segment_align(output_sam_file, fastq_file, fastq_format, seed_length, 
                        multihits, mismatches, num_threads, bowtie_bin, bowtie_index):
    # get the read length to determine how much trimming is needed
    read_length = get_read_length(fastq_file)
    trim3 = read_length - seed_length
    args = [bowtie_bin, "-q", "-S", 
            #"-u", "50000",
            #"--shmem", 
            "-p", str(num_threads),
            "--tryhard",
            "--%s" % fastq_format,
            "-l", str(seed_length),
            "-k", str(multihits),
            "-m", str(multihits),
            "-n", str(mismatches),
            "--trim3", str(trim3)]
    args += [bowtie_index, fastq_file, output_sam_file]
    return args


def align_mate(fastq_file, output_sam_file, seed_length,
               fastq_format, multihits, mismatches, 
               bowtie_threads, bowtie_bin, bowtie_index,
               gene_to_genome_map, tmp_dir):
    try:
        logging.info("Aligning reads from file %s" % (fastq_file))
        # align the mate reads
        aln_fifo_dir, aln_fifo_file = make_fifo(tmp_dir)
        args = setup_segment_align(aln_fifo_file, fastq_file, fastq_format, seed_length, multihits, mismatches, bowtie_threads, bowtie_bin, bowtie_index)
        logging.debug("Bowtie args: %s" % str(args))
        aln_p = subprocess.Popen(args)
        # translate splice reads to genomic reads
        translate_p = multiprocessing.Process(target=transcriptome_to_genome, args=(aln_fifo_file, output_sam_file, gene_to_genome_map))
        translate_p.daemon = True
        translate_p.start()
        # pass to downstream process
        translate_p.join()
        aln_p.wait()
    finally:
        if os.path.exists(aln_fifo_dir):
            os.unlink(aln_fifo_file)
            shutil.rmtree(aln_fifo_dir)

def align(output_sam_file, output_expr_file, 
          fastq_files, seed_length,
          fastq_format, multihits, mismatches, 
          bowtie_threads, bowtie_bin, bowtie_index,
          gene_bed_file, gene_fasta_prefix,
          insert_size_max):
    tmp_dir = os.path.dirname(output_sam_file)
    # build transcriptome -> genome mappings
    logging.info("Constructing table to convert transcriptome references to genome...")
    gene_to_genome_map = build_gene_to_genome_map(gene_bed_file, gene_fasta_prefix)
    try:
        # align mates independently
        processes = []
        fifos = []
        for fastq_file in fastq_files:
            fifo_dir, fifo_file = make_fifo(tmp_dir)
            p = multiprocessing.Process(target=align_mate,
                                        args=(fastq_file, fifo_file, seed_length,
                                              fastq_format, multihits, mismatches, 
                                              bowtie_threads, bowtie_bin, bowtie_index, 
                                              gene_to_genome_map, tmp_dir))
            p.start()
            processes.append(p)
            fifos.append((fifo_dir, fifo_file))
        # merge mate alignments
        sam_fifo_files = [x[1] for x in fifos]
        out_p = multiprocessing.Process(target=join_paired_end_reads, 
                                        args=(sam_fifo_files, output_sam_file, output_expr_file, gene_bed_file, insert_size_max))
        out_p.daemon = True
        out_p.start()
        out_p.join()
        for p in processes:
            p.join()
    finally:
        for fifo_dir, fifo_file in fifos:
            os.unlink(fifo_file)
            shutil.rmtree(fifo_dir)
    return

def check_params(parser, options):
    read_lengths = []
    for fq in options.fastq_files:
        logging.debug("Checking read length for file %s" % (fq))
        read_lengths.append(get_read_length(fq))
        logging.debug("Read length for file %s: %d" % (fq, read_lengths[-1]))
    if any(options.seed_length > rlen for rlen in read_lengths):
        parser.error("seed length %d cannot be longer than read length" % (options.seed_length))

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")    
    parser = argparse.ArgumentParser()
    parser.add_argument("--bowtie-bin", dest="bowtie_bin")
    parser.add_argument("--bowtie-index", dest="bowtie_index")
    parser.add_argument("--bowtie-threads", type=int, dest="bowtie_threads", default=1)
    parser.add_argument("--gene-bed", dest="gene_bed_file")
    parser.add_argument("--gene-fasta-prefix", dest="gene_fasta_prefix")
    parser.add_argument("--multihits", type=int, dest="multihits", default=40)
    parser.add_argument("--mismatches", type=int, dest="mismatches", default=2)
    parser.add_argument("--seed-length", type=int, dest="seed_length", default=25)
    parser.add_argument("--quals", dest="fastq_format")
    parser.add_argument("--insert-size-max", type=int, dest="insert_size_max", default=400)
    parser.add_argument("fastq_files", nargs=2)
    parser.add_argument("output_file")
    parser.add_argument("expr_file")
    options = parser.parse_args()
    check_params(parser, options)        
    align(options.output_file, options.expr_file, 
          options.fastq_files, options.seed_length,
          options.fastq_format, options.multihits, options.mismatches,
          options.bowtie_threads, options.bowtie_bin, options.bowtie_index, 
          options.gene_bed_file, options.gene_fasta_prefix,
          options.insert_size_max)

if __name__ == '__main__': main()