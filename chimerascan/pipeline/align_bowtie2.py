'''
Created on Jun 1, 2012

@author: mkiyer
'''
import os
import sys
import subprocess
import logging

from chimerascan.lib import config
from chimerascan.lib.base import LibraryTypes

import chimerascan.pipeline

_pipeline_dir = chimerascan.pipeline.__path__[0]

_bowtie2_pe_args = ['--phred33', '--end-to-end', '--very-sensitive', 
                    '--reorder', '--no-mixed', '--no-discordant']
_bowtie2_pe_sr_args = ['--phred33', '--end-to-end', '--very-sensitive', 
                       '--reorder']
            
def get_bowtie_library_type(library_type):
    """
    returns the bowtie library type option '--fr' or '--ff' corresponding
    to the first two characters of the library type string
    """
    return library_type[0:2]

def bowtie2_align_transcriptome_pe(transcriptome_index,
                                   genome_index,
                                   transcript_file,                                   
                                   fastq_files,
                                   unaligned_path,
                                   bam_file,
                                   log_file,
                                   library_type,
                                   min_fragment_length=0,
                                   max_fragment_length=1000,
                                   max_transcriptome_hits=1,
                                   num_processors=1):
    """
    align reads to a transcriptome index, convert SAM to BAM,
    and translate alignments to genomic coordinates
    """
    # check num processors
    if num_processors < 2:
        logging.warning("Transcriptome alignment uses a piping approach to "
                        "reduce file I/O and improve overall runtime. This "
                        "requires at least 2 cpu cores. Please set "
                        "num_processors >= 2")
        num_processors = 2
    # setup bowtie2 library type param
    library_type_param = get_bowtie_library_type(library_type)
    # setup bowtie2 command line args
    args = [config.BOWTIE2_BIN]
    args.extend(_bowtie2_pe_args)
    args.extend(['-p', num_processors-1, 
                 '-M', max_transcriptome_hits,
                 '-I', min_fragment_length,
                 '-X', max_fragment_length,
                 '--%s' % library_type_param,
                 '-x', transcriptome_index,
                 '-1', fastq_files[0],
                 '-2', fastq_files[1],
                 '--un-conc', unaligned_path])
    args = map(str, args)
    logging.debug("Alignment args: %s" % (' '.join(args)))    
    # kickoff alignment process
    logfh = open(log_file, "w")
    aln_p = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=logfh)
    # pipe the bowtie SAM output to a transcriptome to genome conversion
    # script that writes a genomic BAM file
    py_script = os.path.join(_pipeline_dir, "transcriptome_to_genome.py")
    args = [sys.executable, py_script, "--library-type", library_type, 
            "--input-sam", "--output-sam", genome_index, transcript_file, 
            "-", "-"]
    args = map(str, args)
    logging.debug("Transcriptome to Genome converter args: %s" % 
                  (' '.join(args)))
    convert_p = subprocess.Popen(args, stdin=aln_p.stdout, stdout=subprocess.PIPE, stderr=logfh)
    # pipe the SAM output to a filter that writes in BAM format
    py_script = os.path.join(_pipeline_dir, "sam_to_bam.py")
    args = [sys.executable, py_script, "-", bam_file] 
    logging.debug("SAM to BAM converter args: %s" % (' '.join(args)))
    sam2bam_p = subprocess.Popen(args, stdin=convert_p.stdout, stderr=logfh)    
    # wait for this to finish
    retcode = sam2bam_p.wait()
    if retcode != 0:
        convert_p.terminate()
        aln_p.terminate()
        return config.JOB_ERROR
    retcode = convert_p.wait()
    if retcode != 0:
        aln_p.terminate()
        return config.JOB_ERROR
    retcode = aln_p.wait()
    if retcode != 0:
        return config.JOB_ERROR        
    logfh.close()
    return retcode

def bowtie2_align_pe(index,
                     fastq_files,
                     unaligned_path,
                     bam_file,
                     log_file,
                     library_type,
                     min_fragment_length=0,
                     max_fragment_length=1000,
                     max_hits=1,
                     num_processors=1):
    # setup bowtie2 library type param
    library_type_param = get_bowtie_library_type(library_type)
    args = [config.BOWTIE2_BIN]
    args.extend(_bowtie2_pe_args)
    args.extend(['-p', num_processors, 
                 '-M', max_hits,
                 '-I', min_fragment_length,
                 '-X', max_fragment_length,
                 '--%s' % library_type_param,
                 '-x', index,
                 '-1', fastq_files[0],
                 '-2', fastq_files[1],
                 '--un-conc', unaligned_path])
    args = map(str, args)
    logging.debug("Alignment args: %s" % (' '.join(args)))    
    # kickoff alignment process
    logfh = open(log_file, "w")
    aln_p = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=logfh)
    # pipe the bowtie SAM output to a filter that writes in BAM format
    py_script = os.path.join(_pipeline_dir, "sam_to_bam.py")
    args = [sys.executable, py_script, "-", bam_file] 
    logging.debug("SAM to BAM converter args: %s" % (' '.join(args)))
    sam2bam_p = subprocess.Popen(args, stdin=aln_p.stdout, stderr=logfh)
    # wait for this to finish
    retcode = sam2bam_p.wait()
    if retcode != 0:
        aln_p.terminate()
        return config.JOB_ERROR
    retcode = aln_p.wait()
    if retcode != 0:
        return config.JOB_ERROR        
    logfh.close()
    return retcode

def parse_fastq(line_iter):
    with line_iter:
        while True:
            lines = [line_iter.next().rstrip() for x in xrange(4)]
            yield lines

def trim_and_merge_fastq(infiles, outfh, trimmed_outfh, segment_length):
    fqiters = [parse_fastq(open(f)) for f in infiles]    
    try:
        while True:
            pe_lines = [fqiter.next() for fqiter in fqiters]
            for readnum,lines in enumerate(pe_lines):
                # encode a '0' or '1' as the first character of the line
                lines[0] = "@%d%s" % (readnum, lines[0][1:])
                seqlen = len(lines[1])
                # output full length reads
                print >>outfh, '\n'.join(lines)
                # output trimmed reads                
                if seqlen > segment_length:
                    lines[1] = lines[1][:segment_length]
                    lines[3] = lines[3][:segment_length]
                print >>trimmed_outfh, '\n'.join(lines)
    except StopIteration:
        pass
    return config.JOB_SUCCESS

def bowtie2_align_pe_sr(index,
                        transcript_file,
                        fastq_files,
                        bam_file,
                        log_file,
                        tmp_dir,
                        segment_length,
                        max_hits=1,
                        num_processors=1):
    # create tmp interleaved fastq file
    interleaved_fastq_file = os.path.join(tmp_dir, config.INTERLEAVED_FASTQ_FILE)
    interleaved_trimmed_fastq_file = os.path.join(tmp_dir, config.INTERLEAVED_TRIMMED_FASTQ_FILE)
    fqfh = open(interleaved_fastq_file, 'w')
    trimmed_fqfh = open(interleaved_trimmed_fastq_file, 'w')
    retcode = trim_and_merge_fastq(fastq_files, 
                                   fqfh, trimmed_fqfh,
                                   segment_length)
    fqfh.close()
    trimmed_fqfh.close()
    if retcode != config.JOB_SUCCESS:
        logging.error("[FAILED] trimming FASTQ files")
        if os.path.exists(interleaved_fastq_file):
            os.remove(interleaved_fastq_file)
        if os.path.exists(interleaved_trimmed_fastq_file):
            os.remove(interleaved_trimmed_fastq_file)
        return config.JOB_ERROR
    #
    # Align trimmed reads
    #
    logfh = open(log_file, "w")
    # setup bowtie2 library type param
    args = [config.BOWTIE2_BIN]
    args.extend(_bowtie2_pe_sr_args)
    args.extend(['-p', num_processors, 
                 '-k', max_hits,
                 '-x', index,
                 '-U', interleaved_trimmed_fastq_file])
    args = map(str, args)
    logging.debug("Alignment args: %s" % (' '.join(args)))
    aln_p = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=logfh)
    #
    # Convert to BAM, also extend sequences back to full length by 
    # adding padding to CIGAR string
    #
    py_script = os.path.join(_pipeline_dir, "sam_to_bam_pesr.py")
    args = [sys.executable, py_script, "-", interleaved_fastq_file, bam_file]
    logging.debug("SAM to BAM converter args: %s" % (' '.join(args)))
    sam2bam_p = subprocess.Popen(args, stdin=aln_p.stdout, stderr=logfh)
    retcode = sam2bam_p.wait()
    if retcode != 0:
        logging.debug("Error during SAM to BAM conversion")
        aln_p.terminate()
        return config.JOB_ERROR
    retcode = aln_p.wait()
    if retcode != 0:
        logging.debug("Error during alignment")
        return config.JOB_ERROR        
    logfh.close()
    return retcode

def bowtie2_align_local(transcriptome_index,
                        genome_index,
                        transcript_file,                                   
                        fastq_file,
                        bam_file,
                        log_file,
                        local_anchor_length,
                        local_multihits,
                        num_processors=1):
    """
    align reads to a transcriptome index, convert SAM to BAM,
    and translate alignments to genomic coordinates
    """
    # check num processors
    if num_processors < 2:
        logging.warning("Transcriptome alignment uses a piping approach to "
                        "reduce file I/O and improve overall runtime. This "
                        "requires at least 2 cpu cores. Please set "
                        "num_processors >= 2")
        num_processors = 2
    # seed substrings should be size of local anchor length
    seed_length = local_anchor_length
    # minimum score to report and alignment is twice local anchor length
    score_min = 2*local_anchor_length
    # setup bowtie2 command line args
    args = [config.BOWTIE2_BIN,
            '-q',
            '--phred33',
            '-N', '0',
            '-L', seed_length,
            '-i', 'C,1,0',
            '--local',
            '--score-min', 'C,%d,0' % (score_min),
            '-k', local_multihits,
            '-R', '2',
            '-p', num_processors-1,             
            '--reorder',
            '-x', transcriptome_index,
            '-U', fastq_file]
    args = map(str, args)
    logging.debug("Alignment args: %s" % (' '.join(args)))    
    # kickoff alignment process
    logfh = open(log_file, "w")
    aln_p = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=logfh)
    # pipe the bowtie SAM output to a transcriptome to genome conversion
    # script that writes a genomic BAM file
    py_script = os.path.join(_pipeline_dir, "transcriptome_to_genome.py")
    args = [sys.executable, py_script, 
            "--library-type", LibraryTypes.FR_UNSTRANDED,
            "--input-sam", "--output-sam", genome_index, transcript_file, 
            "-", "-"]
    args = map(str, args)
    logging.debug("Transcriptome to Genome converter args: %s" % 
                  (' '.join(args)))
    convert_p = subprocess.Popen(args, stdin=aln_p.stdout, stdout=subprocess.PIPE, stderr=logfh)
    # pipe the SAM output to a filter that writes in BAM format
    py_script = os.path.join(_pipeline_dir, "sam_to_bam.py")
    args = [sys.executable, py_script, "-", bam_file] 
    logging.debug("SAM to BAM converter args: %s" % (' '.join(args)))
    sam2bam_p = subprocess.Popen(args, stdin=convert_p.stdout, stderr=logfh)    
    # wait for this to finish
    retcode = sam2bam_p.wait()
    if retcode != 0:
        convert_p.terminate()
        aln_p.terminate()
        return config.JOB_ERROR
    retcode = convert_p.wait()
    if retcode != 0:
        aln_p.terminate()
        return config.JOB_ERROR
    retcode = aln_p.wait()
    if retcode != 0:
        return config.JOB_ERROR        
    logfh.close()
    return retcode


def bowtie2_align_sr(index,
                     fastq_file,
                     bam_file,
                     log_file,
                     library_type,
                     maxhits=1,
                     num_processors=1):
    # align reads
    logfh = open(log_file, "w")
    # setup bowtie2 library type param
    args = [config.BOWTIE2_BIN, 
            '-p', num_processors, '--phred33',
            '--end-to-end', '--very-sensitive', '--reorder',
            '-M', maxhits,
            '-x', index,
            '-U', fastq_file]
    args = map(str, args)
    logging.debug("Alignment args: %s" % (' '.join(args)))
    aln_p = subprocess.Popen(args, stdout=subprocess.PIPE, stderr=logfh)
    # pipe the bowtie SAM output to a filter that writes in BAM format
    py_script = os.path.join(_pipeline_dir, "sam_to_bam.py")
    args = [sys.executable, py_script, "-", bam_file] 
    logging.debug("SAM to BAM converter args: %s" % (' '.join(args)))
    sam2bam_p = subprocess.Popen(args, stdin=aln_p.stdout, stderr=logfh)
    # wait for this to finish
    retcode = sam2bam_p.wait()
    if retcode != 0:
        if os.path.exists(bam_file):
            os.remove(bam_file)
        aln_p.terminate()
    else:
        retcode = aln_p.wait()
        if retcode != 0:
            if os.path.exists(bam_file):
                os.remove(bam_file)
    logfh.close()
    return retcode