'''
Created on Jan 7, 2011

@author: mkiyer
'''
'''
Created on Oct 22, 2010

@author: mkiyer
'''
import sys
import os
import logging
import subprocess

from config import MIN_SEGMENT_LENGTH, JOB_SUCCESS
from base import get_read_length

class AlignError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

def determine_read_segments(read_length, segment_length, segment_trim, trim5, trim3):
    # figure out how many segments are available
    trimmed_rlen = read_length - trim5 - trim3
    num_segments = trimmed_rlen / segment_length    
    if segment_trim:
        read_end = segment_length * num_segments
    else:
        read_end = read_length - trim3    
    if num_segments == 1:
        # if just one segment is available, just use a single segment length
        return [(trim5, read_end)]
    start = trim5 + segment_length
    segments = [(trim5, start)]
    for x in xrange(1, num_segments - 1):
        segments.append((start, start + segment_length))
        start += segment_length
    if start < read_end:        
        segments.append((start, read_end))
    return segments

def align_segments(fastq_files, output_bam_file, segments,
                   fastq_format, multihits, mismatches, 
                   num_threads, bowtie_bin, bowtie_index,
                   bowtie_mode):    
    #
    # Cut reads into segments and merge paired-end reads 
    # into a single FASTQ file
    #
    seg_starts_str = ','.join(map(str, [s[0] for s in segments]))
    seg_ends_str = ','.join(map(str, [s[1] for s in segments]))    
    py_script = os.path.join(os.path.dirname(__file__), "segment_reads.py")
    args = [sys.executable, py_script, seg_starts_str, seg_ends_str]
    args.extend(fastq_files)
    logging.debug("Read segmentation args: %s" % (' '.join(args)))
    seg_p = subprocess.Popen(args, stdout=subprocess.PIPE)
    max_segment_length = max((s[1] - s[0]) for s in segments)
    #
    # Align the segmented reads
    #
    args = [bowtie_bin, "-q", "-S", 
            "-p", str(num_threads),
            "--tryhard",
            "--%s" % fastq_format,
            "-k", str(multihits),
            "-m", str(multihits),
            bowtie_mode, str(mismatches)]
    if bowtie_mode == "-n":
        args.extend(["-l", str(max_segment_length)])
    #args += [bowtie_index, "-", output_sam_file]
    args += [bowtie_index, "-"]
    logging.debug("Alignment args: %s" % (' '.join(args)))
    aln_p = subprocess.Popen(args, stdin=seg_p.stdout, stdout=subprocess.PIPE)
    #
    # Merge segmented alignments
    #
    py_script = os.path.join(os.path.dirname(__file__), "join_segmented_alignments.py")
    args = [sys.executable, py_script]
    if len(fastq_files) == 1:
        args.append("--sr")
    args.extend(["-", fastq_files[0], output_bam_file])
    logging.debug("Join segmented alignments args: %s" % (' '.join(args)))
    out_p = subprocess.Popen(args, stdin=aln_p.stdout)
    out_p.wait()    
    aln_p.wait()
    seg_p.wait()

def check_fastq_files(fastq_files, segment_length, trim5, trim3):
    # check that input fastq files exist
    read_lengths = []
    for mate,fastq_file in enumerate(fastq_files):
        if not os.path.isfile(fastq_file):
            raise AlignError("mate '%d' fastq file '%s' is not valid" % 
                             (mate, fastq_file))
        logging.debug("Checking read length for file %s" % (fastq_file))
        read_lengths.append(get_read_length(fastq_file))
        logging.debug("Read length for file %s: %d" % 
                      (fastq_file, read_lengths[-1]))
    # check that mate read lengths are equal
    if len(set(read_lengths)) > 1:
        logging.error("read lengths mate1=%d and mate2=%d are unequal" % 
                      (read_lengths[0], read_lengths[1]))
        return False
    rlen = read_lengths[0]
    trimmed_rlen = rlen - trim5 - trim3
    # check that segment length >= MIN_SEGMENT_LENGTH 
    if segment_length < MIN_SEGMENT_LENGTH:
        raise AlignError("segment length (%d) too small (min is %d)" % 
                         (segment_length, MIN_SEGMENT_LENGTH))
    # check that segment length < trimmed read length
    if segment_length > trimmed_rlen:
        raise AlignError("segment length (%d) longer than trimmed read length (%d)" % 
                         (segment_length, trimmed_rlen))
    return read_lengths[0]

def align(fastq_files, fastq_format, 
          bowtie_index, output_bam_file,  
          bowtie_bin="bowtie",
          num_processors=2,
          segment_length=25, 
          segment_trim=False,
          trim5=0, 
          trim3=0, 
          multihits=40, 
          mismatches=2,
          bowtie_mode="-n"):
    # check fastq files
    read_length = check_fastq_files(fastq_files, segment_length, trim5, trim3)
    # divide reads into segments
    # TODO: trim reads with polyA tails or bad quality scores
    # signify mate1 and mate2
    segments = determine_read_segments(read_length, segment_length, segment_trim, trim5, trim3)
    logging.info("Dividing %dbp reads into %d segments: %s" %
                 (read_length, len(segments), segments))      
    # run paired-end segmented aligner
    logging.info("Running segmented alignment")
    align_segments(fastq_files, output_bam_file, segments,
                   fastq_format, multihits, mismatches, num_processors,
                   bowtie_bin, bowtie_index, bowtie_mode)
    logging.info("Alignment completed")
    return JOB_SUCCESS

def main():
    from optparse import OptionParser
    import config
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = OptionParser("usage: %prog [options] <output_dir> <mate1.fq> [<mate2.fq>]")
    parser.add_option("-p", "--processors", dest="num_processors", 
                      type="int", default=1)
    parser.add_option("--bowtie-bin", dest="bowtie_bin", default="bowtie", 
                      help="Path to 'bowtie' program")
    parser.add_option("--bowtie-index", dest="bowtie_index",
                      help="Path to bowtie index")
    parser.add_option("--bowtie-mode-v", action="store_true", 
                      dest="bowtie_mode_v", default=False,
                      help="Run bowtie with -v to ignore quality scores")
    parser.add_option("--multihits", type="int", dest="multihits", 
                      default=100)
    parser.add_option("--mismatches", type="int", dest="mismatches", 
                      default=2)
    parser.add_option("--segment-length", type="int", dest="segment_length", 
                      default=25)
    parser.add_option("--trim5", type="int", dest="trim5", 
                      default=0)
    parser.add_option("--trim3", type="int", dest="trim3", 
                      default=0)
    parser.add_option("--segment-trim", action="store_true", default=False,
                      help="Trim reads to be an exact multiple of the " 
                      "segment length")
    parser.add_option("--quals", dest="fastq_format", default="phred33-quals")
    options, args = parser.parse_args()
    # extract command line arguments
    output_bam_file = args[0]
    fastq_files = args[1:]
    if options.bowtie_mode_v:
        bowtie_mode = "-v"
    else:
        bowtie_mode = "-n" 
    retcode = align(fastq_files, options.fastq_format,
                    options.bowtie_index, 
                    output_bam_file,  
                    bowtie_bin=options.bowtie_bin,
                    num_processors=options.num_processors,
                    segment_length=options.segment_length,
                    segment_trim=options.segment_trim,
                    trim5=options.trim5,
                    trim3=options.trim3,
                    multihits=options.multihits,
                    mismatches=options.mismatches,
                    bowtie_mode=bowtie_mode)
    sys.exit(retcode) 

if __name__ == '__main__':
    main()
