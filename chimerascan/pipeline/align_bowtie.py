'''
Created on Jan 22, 2011

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
import logging
import subprocess
import sys

# local imports
import chimerascan.pysam as pysam
from chimerascan.lib.base import get_read_length
from fix_alignment_ordering import fix_pe_alignment_ordering, fix_sr_alignment_ordering

translate_quals = {'solexa': '--solexa-quals',
                   'illumina': '--solexa1.3-quals',
                   'sanger': '--phred33-quals'}

def align_pe_full(fastq_files, 
                  bowtie_index,
                  output_bam_file, 
                  unaligned_fastq_param=None,
                  maxmultimap_fastq_param=None,
                  min_fragment_length=0,
                  max_fragment_length=1000,
                  trim5=0,
                  trim3=0,
                  library_type="fr",
                  num_processors=1, 
                  quals="sanger", 
                  multihits=100, 
                  mismatches=2, 
                  bowtie_bin="bowtie", 
                  bowtie_mode="-n",
                  log_file=None):
    # convert from quality format description strings to
    # options used to run bowtie
    fastq_format_option = translate_quals[quals]
    # the first two characters of the library type option
    # should be either "fr", "rf", or "ff" in order to match
    # the bowtie orientation options
    read_orientation = library_type[:2]        
    args = [bowtie_bin, "-q", "-S", 
            "-p", str(num_processors),
            fastq_format_option,
            "-k", str(multihits),
            "-m", str(multihits),
            bowtie_mode, str(mismatches),
            "--minins", min_fragment_length,
            "--maxins", max_fragment_length,
            "--trim5", trim5,
            "--trim3", trim3,
            "--%s" % read_orientation]
    if unaligned_fastq_param is not None:
        args.extend(["--un", unaligned_fastq_param])
    if maxmultimap_fastq_param is not None:
        args.extend(["--max", maxmultimap_fastq_param])
    # use the entire read length as the "seed" here
    if bowtie_mode == "-n":
        read_length = get_read_length(fastq_files[0])
        args.extend(["-l", str(read_length)])
    args += [bowtie_index, 
             "-1", fastq_files[0],
             "-2", fastq_files[1]]
    #aligned_sam_file]
    args = map(str, args)
    logging.debug("Bowtie alignment args: %s" % (' '.join(args)))
    aln_p = subprocess.Popen(args, stdout=subprocess.PIPE)
    # pipe the bowtie SAM output to a filter that writes BAM format
    args = [sys.executable, __file__, "--multihits", str(multihits),
            output_bam_file, fastq_files[0]]
    logging.debug("SAM to BAM converter args: %s" % (' '.join(args)))
    if log_file is not None:
        logfh = open(log_file, "w")
    else:
        logfh = None    
    retcode = subprocess.call(args, stdin=aln_p.stdout, stderr=logfh)
    if logfh is not None:
        logfh.close()
    if retcode != 0:
        return retcode
    return aln_p.wait()

def align_sr_full(fastq_file, 
                  bowtie_index,
                  output_bam_file, 
                  trim5=0,
                  trim3=0,
                  num_processors=1, 
                  quals="sanger",
                  multihits=100, 
                  mismatches=2, 
                  bowtie_bin="bowtie", 
                  bowtie_mode="-n",
                  log_file=None):               
    read_length = get_read_length(fastq_file)
    fastq_format_option = translate_quals[quals]
    args = [bowtie_bin, "-q", "-S", 
            "-p", str(num_processors),
            fastq_format_option,
            "-k", str(multihits),
            "-m", str(multihits),
            bowtie_mode, str(mismatches),
            "--trim5", trim5,
            "--trim3", trim3]
    # use the entire read length as the "seed" here
    if bowtie_mode == "-n":
        args.extend(["-l", str(read_length)])
    args += [bowtie_index, fastq_file]
    #aligned_sam_file]
    args = map(str, args)
    logging.debug("Bowtie alignment args: %s" % (' '.join(args)))
    aln_p = subprocess.Popen(args, stdout=subprocess.PIPE)
    # pipe the bowtie SAM output to a filter that writes BAM format
    args = [sys.executable, __file__, "--multihits", str(multihits), "--sr",
            output_bam_file, fastq_file]
    logging.debug("SAM to BAM converter args: %s" % (' '.join(args)))
    if log_file is not None:
        logfh = open(log_file, "w")
    else:
        logfh = None
    retcode = subprocess.call(args, stdin=aln_p.stdout, stderr=logfh)
    if logfh is not None:
        logfh.close()
    if retcode != 0:
        return retcode
    return aln_p.wait()


def sam_stdin_to_bam(output_bam_file, input_fastq_file, multihits, 
                     is_paired=True, keep_unmapped=True):
    samfh = pysam.Samfile("-", "r")
    bamfh = pysam.Samfile(output_bam_file, "wb", template=samfh)    
    num_unmapped = 0
    num_multihits = 0
    if is_paired:
        for pe_reads in fix_pe_alignment_ordering(samfh, 
                                                  open(input_fastq_file), 
                                                  is_paired=is_paired):
            for reads in pe_reads:
                for r in reads:
                    if r.is_unmapped:
                        xm_tag = r.opt('XM')
                        if xm_tag < multihits:
                            num_unmapped += 1
                            if not keep_unmapped:
                                continue
                        num_multihits += 1
                    bamfh.write(r)
    else:
        for reads in fix_sr_alignment_ordering(samfh, open(input_fastq_file)): 
            for r in reads:
                if r.is_unmapped:
                    xm_tag = r.opt('XM')
                    if xm_tag < multihits:
                        num_unmapped += 1
                        if not keep_unmapped:
                            continue
                    num_multihits += 1
                bamfh.write(r)
    bamfh.close()
    samfh.close()
    logging.debug("[SAMTOBAM] Filtered %d unmapped reads" % (num_unmapped))
    logging.debug("[SAMTOBAM] Found %d multimapping (>%d) reads" % 
                  (num_multihits, multihits))
    logging.info("[SAMTOBAM] Finished converting SAM -> BAM")

if __name__ == '__main__':
    from optparse import OptionParser
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = OptionParser("usage: %prog [options] <chimeras.bedpe> <out.fasta> <out.juncs>")
    parser.add_option("--multihits", type="int", dest="multihits", default=100)
    parser.add_option("--sr", action="store_false", dest="is_paired", default=True)
    parser.add_option("--un", action="store_true", dest="keep_unmapped", default=False)
    options, args = parser.parse_args()
    output_bam_file = args[0]
    input_fastq_file = args[1]
    sam_stdin_to_bam(output_bam_file, input_fastq_file, 
                     multihits=options.multihits,
                     is_paired=options.is_paired,
                     keep_unmapped=options.keep_unmapped)
