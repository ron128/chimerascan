'''
Created on Oct 14, 2010

@author: mkiyer
'''
import argparse
import tempfile
import os
import sys
import logging
import collections
import subprocess
import multiprocessing
import shutil
import operator
import pysam

from base import get_read_length

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

def setup_bowtie_segment_align(output_sam_file, fastq_file, fastq_format, 
                               seed_length, multihits, mismatches, num_threads, 
                               bowtie_bin, bowtie_index):
    args = [bowtie_bin, "-q", "-S", 
#            "--shmem",
            "-p", str(num_threads),
            "--tryhard",
            "--%s" % fastq_format,
            "-l", str(seed_length),
            "-k", str(multihits),
            "-m", str(multihits),
            "-n", str(mismatches)]
    args += [bowtie_index, fastq_file, output_sam_file]
    return args

def align_segment(fastq_file, output_sam_file, segment_start, segment_end,
                  fastq_format, multihits, mismatches, 
                  num_threads, bowtie_bin, bowtie_index,
                  tmp_dir):    
    # cut read into segments
    py_script = os.path.join(os.path.dirname(__file__), "segment_reads.py")
    args = [sys.executable, py_script, fastq_file, str(segment_start), str(segment_end)]
    seg_p = subprocess.Popen(args, stdout=subprocess.PIPE)
    # align the segmented reads
    segment_length = segment_end - segment_start
    args = setup_bowtie_segment_align(output_sam_file, "-", fastq_format, 
                                      segment_length, multihits, mismatches, num_threads, 
                                      bowtie_bin, bowtie_index)
    aln_p = subprocess.Popen(args, stdin=seg_p.stdout)
    aln_p.wait()
    seg_p.wait()

def iter_multihit_alignments(samfh):    
    reads = []
    for read in samfh:
        if len(reads) > 0 and read.qname != reads[-1].qname:
            yield reads
            reads = []
        reads.append(read)
    if len(reads) > 0:
        yield reads

def parse_segment_sam_files(samfhs, maxlen=100000):
    buf = []
    for i in xrange(maxlen):
        buf.append([list() for i in xrange(len(samfhs))])
    buf_size = 0
    next_buf_ind = 0
    buf_ind = 0
    qname_ind_map = {}
    sam_iters = [iter_multihit_alignments(f) for f in samfhs]
    #for read_segments in itertools.izip(*sam_iters):
    try:
        while True:
            read_segments = [sam_iter.next() for sam_iter in sam_iters]
            for segment_num,reads in enumerate(read_segments):
                for read in reads:
                    # set read flags
                    read.is_read1 = True
                    read.qname = read.qname[:read.qname.rfind("_")]
                    read.is_paired = False
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
                        buf[buf_ind] = [list() for i in xrange(len(samfhs))]
                    else:
                        buf_ind = qname_ind_map[read.qname]
                    buf[buf_ind][segment_num].append(read)
    except StopIteration:
        pass
    for buf_ind in xrange(buf_size):
        yield buf[buf_ind]

def build_segment_alignment_dict(aln_dict, reads, seg_ind, segment, junc_pos_map):
    for i,r in enumerate(reads):
        if r.is_unmapped:
            continue
        junc_pos = junc_pos_map[r.rname]
        is_spanning = (r.pos < junc_pos) and (r.aend > junc_pos)
        offset = r.pos + segment[0] if r.is_reverse else r.pos - segment[0]
        aln_dict[(r.rname, r.is_reverse, offset)].append((seg_ind, i, is_spanning))

def get_consecutive_spanning_segments(hits, segment_alignments):
    begin_ind = None
    end_ind = None
    segs_spanning = False
    segs = []
    for seg_ind, read_ind, is_spanning in sorted(hits, key=operator.itemgetter(0)):    
        if begin_ind is None:
            begin_ind = seg_ind
        elif seg_ind != end_ind:
            if (end_ind - begin_ind >= 2) and segs_spanning:
                yield segs 
            begin_ind = seg_ind
            segs_spanning = False
            segs = []            
        end_ind = seg_ind + 1
        segs.append((seg_ind, segment_alignments[seg_ind][read_ind]))
        segs_spanning = segs_spanning or is_spanning
    if (end_ind - begin_ind >= 2) and segs_spanning:
        yield segs

def find_valid_segment_alignments(segment_alignments, segments, junc_pos_map):
    aln_dict = collections.defaultdict(lambda: [])
    for i in xrange(len(segments)):
        build_segment_alignment_dict(aln_dict, segment_alignments[i], i, segments[i], junc_pos_map)
    for key,hits in aln_dict.iteritems():        
        hit_reads = collections.defaultdict(lambda: [])
        for segs in get_consecutive_spanning_segments(hits, segment_alignments):
            hit_reads[len(segs)].append(segs) 
        # if there are multiple hits where different numbers of segments
        # map, then just take the hits where the most number of segments
        # successfully map
        if len(hit_reads) > 0:        
            best_num_segs = sorted(hit_reads.iterkeys(), reverse=True)[0]
            for segs in hit_reads[best_num_segs]:
                yield segs

def parse_MD_tag(val):
    x = 0
    mdops = []
    for y in xrange(len(val)):
        if val[y].isalpha():
            offset = int(val[x:y])
            base = val[y]
            mdops.append((offset, base))
            x = y + 1
#    if x < len(val):
#        mdops.append(int(val[x:]))
    return mdops

def create_bowtie_output(reads, segments, references):
    seqlen = segments[-1][1] - segments[0][0]    
    seq = ['N'] * seqlen
    qual = ['!'] * seqlen
    strand = "-" if reads[0].is_reverse else "+"
    pos = reads[0].pos if strand == "+" else reads[-1].pos
    mismatches = set()
    for read in reads:
        newpos = read.pos - pos
        newend = read.aend - pos
        seq[newpos:newend] = read.seq
        qual[newpos:newend] = read.qual
        # parse mismatches
        current_offset = 0
        for offset,base in parse_MD_tag(read.opt("MD")):
            current_offset += offset
            mdpos = newpos + current_offset
            current_offset += 1
            mismatches.add((mdpos, base, seq[mdpos]))   
    mismatch_string = ','.join("%d:%s>%s" % (mdpos, refbase, base) for mdpos, refbase, base in sorted(mismatches, key=operator.itemgetter(0))) 
    fields = [reads[0].qname, strand,
              references[reads[0].rname],
              str(pos),
              ''.join(seq),
              ''.join(qual),
              "0",
              mismatch_string]
    return '\t'.join(fields)

def join_segments(segment_sam_files, segments, output_file):
    samfhs = [pysam.Samfile(f, "r") for f in segment_sam_files]
    references = samfhs[0].references
    # build map of where the fusion junction occurs on
    # each reference sequence
    junc_pos_map = {}
    for i,rname in enumerate(references):
        left_length = int(rname.split(":")[1])
        junc_pos_map[i] = left_length
    # find combinations of junction spanning
    outfh = open(output_file, "w")
    for segment_alignments in parse_segment_sam_files(samfhs):        
        for segs in find_valid_segment_alignments(segment_alignments, segments, junc_pos_map):
            reads = []
            intervals = []
            # get offset distance from original read
            offset = segments[segs[0][0]][0]
            for ind,r in segs:
                reads.append(r)                
                intervals.append((segments[ind][0] - offset, segments[ind][1] - offset))
            print >>outfh, create_bowtie_output(reads, intervals, references)
    outfh.close()

def determine_read_segments(read_length, segment_length):
    num_segments = int(1 + (read_length / float(segment_length)))
    segment_offset = (read_length - segment_length) / float(num_segments - 1)
    segments = [(0, segment_length)]
    for i in xrange(1, num_segments-1):
        start = int(i * segment_offset)        
        segments.append((start, start + segment_length))
    segments.append((read_length - segment_length, read_length))
    return segments

def segmented_aligner(output_file, fastq_file, fastq_format, 
                      segment_length, multihits, mismatches, 
                      bowtie_threads, bowtie_bin, bowtie_index):
    tmp_dir = os.path.dirname(output_file)
    # determine how to break read into segments
    read_length = get_read_length(fastq_file)
    segments = determine_read_segments(read_length, segment_length)
    logging.info("Dividing %dbp reads into %d segments: %s" %
                 (read_length, len(segments), segments))    
    try:
        # align segments independently 
        seg_processes = []
        seg_fifos = []
        for segment in segments:        
            fifo_dir, fifo_file = make_fifo(tmp_dir)
            p = multiprocessing.Process(target=align_segment, 
                                        args=(fastq_file, fifo_file, segment[0], segment[1],
                                              fastq_format, multihits, mismatches, 
                                              bowtie_threads, bowtie_bin, bowtie_index, tmp_dir))
            p.start()
            seg_processes.append(p)
            seg_fifos.append((fifo_dir, fifo_file))
        # merge segment alignments
        sam_fifo_files = [x[1] for x in seg_fifos]
        out_p = multiprocessing.Process(target=join_segments, args=(sam_fifo_files, segments, output_file))
        out_p.daemon = True
        out_p.start()
        out_p.join()
        for p in seg_processes:
            p.join()
    finally:
        for fifo_dir, fifo_file in seg_fifos:
            os.unlink(fifo_file)
            shutil.rmtree(fifo_dir)
    return

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")    
    parser = argparse.ArgumentParser()
    parser.add_argument("--bowtie-bin", dest="bowtie_bin")
    parser.add_argument("--bowtie-index", dest="bowtie_index")
    parser.add_argument("--bowtie-threads", dest="bowtie_threads", default=1)
    parser.add_argument("--segment-multihits", type=int, dest="multihits", default=40)
    parser.add_argument("--segment-mismatches", type=int, dest="mismatches", default=2)
    parser.add_argument("--segment-length", type=int, dest="segment_length", default=-1)
    parser.add_argument("--quals", dest="fastq_format")
    parser.add_argument("fastq_file")
    parser.add_argument("output_file")
    options = parser.parse_args()
    segmented_aligner(options.output_file, options.fastq_file, options.fastq_format,
                      options.segment_length, options.multihits, options.mismatches,
                      options.bowtie_threads, options.bowtie_bin, options.bowtie_index)

if __name__ == '__main__': main()