'''
Created on Jun 6, 2012

@author: mkiyer
'''
import sys
import os
import logging
import argparse
import collections
import subprocess

from chimerascan import pysam
from chimerascan.lib import config
from chimerascan.lib.base import check_executable, LibraryTypes
from chimerascan.lib.feature import TranscriptFeature
from chimerascan.lib.sam import copy_read, parse_pe_reads, group_read_pairs, pair_reads

def get_sam_header_from_bowtie2_index(index):
    # extract sequence names and lengths from bowtie reference
    args = [config.BOWTIE2_INSPECT_BIN, "-s", index]
    p = subprocess.Popen(args, stdout=subprocess.PIPE)
    stdoutdata = p.communicate()[0]
    retcode = p.wait()
    if retcode != 0:
        return config.JOB_ERROR, {}
    # parse index information and build SAM header
    headerdict = {}
    sqlist = []
    for line in stdoutdata.split('\n'):
        if not line:
            continue
        line = line.strip()
        if not line:
            continue
        if not line.startswith("Sequence"):
            continue
        fields = line.strip().split('\t')
        seqname = fields[1]
        seqlen = int(fields[2])
        sqlist.append({'LN': seqlen, 'SN': seqname})
    headerdict['SQ'] = sqlist
    return headerdict

#
# constants used for CIGAR alignments
#
CIGAR_M = 0 #match  Alignment match (can be a sequence match or mismatch)
CIGAR_I = 1 #insertion  Insertion to the reference
CIGAR_D = 2 #deletion  Deletion from the reference
CIGAR_N = 3 #skip  Skipped region from the reference
CIGAR_S = 4 #softclip  Soft clip on the read (clipped sequence present in <seq>)
CIGAR_H = 5 #hardclip  Hard clip on the read (clipped sequence NOT present in <seq>)
CIGAR_P = 6 #padding  Padding (silent deletion from the padded reference sequence)
CIGAR_E = 7 # sequence match
CIGAR_X = 8  # sequence mismatch
REF_ADVANCING_CIGAR_CODES = frozenset((CIGAR_M, CIGAR_D, CIGAR_N, CIGAR_E, CIGAR_X))

def convert_pos(pos, negstrand, exons):
    """
    find position along transcript
    """
    eindex = 0
    testart = 0     
    exon_size = exons[eindex][1] - exons[eindex][0]
    while pos > (testart + exon_size):
        # skip to next exons
        testart += exon_size
        eindex += 1
        exon_size = exons[eindex][1] - exons[eindex][0]
    toffset = (pos - testart)
    if negstrand:
        newpos = exons[eindex][0] + exon_size - toffset - 1
    else:
        newpos = exons[eindex][0] + toffset
    return newpos, eindex, testart, toffset

def convert_cigar(cigar, negstrand, exons, eindex, testart, toffset):
    exon_size = exons[eindex][1] - exons[eindex][0]    
    newcigar = []
    total_cigar_bp = 0
    for cigarcode, cigarbp in cigar:
        if cigarcode in REF_ADVANCING_CIGAR_CODES:
            # process the aligned cigar bp
            while cigarbp > (exon_size - toffset):
                # subtract remainder of exon from cigar bp
                cigarbp -= (exon_size - toffset)
                # add cigar for remainder of this exon
                newcigar.append((cigarcode, exon_size - toffset))
                total_cigar_bp += (exon_size - toffset)                
                # insert skip CIGAR operation for exon junction
                prev_estart, prev_eend = exons[eindex]
                eindex += 1
                if negstrand:
                    intron_size = (prev_estart - exons[eindex][1])
                else:
                    intron_size = (exons[eindex][0] - prev_eend)
                newcigar.append((CIGAR_N, intron_size))
                total_cigar_bp += intron_size                
                # advance to next exon
                testart += exon_size
                toffset = 0
                exon_size = exons[eindex][1] - exons[eindex][0]     
            # update offset into current exon
            toffset += cigarbp
        # can finish remaining cigarbp within current exon
        newcigar.append((cigarcode, cigarbp))
        total_cigar_bp += cigarbp
    if negstrand:
        # flip CIGAR list
        newcigar.reverse() 
    return newcigar, total_cigar_bp

def convert_read(r, tid_tx_map, genome_rname_tid_map):
    # copy original read and modify it
    newr = copy_read(r)
    if not r.is_unmapped:
        chrom, negstrand, exons = tid_tx_map[r.rname]
        # TODO: remove assert statement
        assert chrom in genome_rname_tid_map
        # convert reference indexes from transcriptome to genome
        genome_tid = genome_rname_tid_map[chrom]
        # find genomic start position of transcript
        newpos, eindex, testart, toffset = convert_pos(r.pos, negstrand, exons)
        # parse and convert transcript cigar string
        newcigar, total_cigar_bp = convert_cigar(r.cigar, negstrand, exons, \
                                                 eindex, testart, toffset)            
        newr.rname = genome_tid
        newr.cigar = newcigar
        # find genomic position of alignment
        if negstrand:
            # set position to left end of transcript
            newr.pos = newpos - total_cigar_bp            
            # flip is_reverse flag
            newr.is_reverse = (not newr.is_reverse)
        else:
            newr.pos = newpos
    return newr

def _setup_and_open_files(genome_index, transcript_feature_file,
                          input_file, output_bam_file, library_type,
                          input_is_sam):
    # create SAM header from genome index
    if not check_executable(config.BOWTIE2_INSPECT_BIN):
        logging.error("Cannot find bowtie2-inspect binary")
        return config.JOB_ERROR
    headerdict = get_sam_header_from_bowtie2_index(genome_index)
    # open input BAM file and add to header
    if input_is_sam:
        mode = "r"
    else:
        mode = "rb"
    infh = pysam.Samfile(input_file, mode)
    inheader = infh.header
    headerdict['HD'] = inheader['HD']
    headerdict['PG'] = inheader['PG']
    # open output BAM file with new header 
    outbamfh = pysam.Samfile(output_bam_file, "wb", header=headerdict)
    # setup reference name mappings
    genome_rname_tid_map = dict((rname,i) for i,rname in enumerate(outbamfh.references))
    tx_rname_tid_map = dict((rname,i) for i,rname in enumerate(infh.references))
    # read transcript feature and prepare data structure for conversion
    logging.debug("Reading transcript features")
    tid_tx_map = {}
    for t in TranscriptFeature.parse(open(transcript_feature_file)):
        exons = [(start, end) for start, end in t.exons]
        negstrand = 1 if t.strand == "-" else 0
        if negstrand:
            exons.reverse()
        tid = tx_rname_tid_map[str(t.tx_id)]
        tid_tx_map[tid] = (t.chrom, negstrand, exons)        
    return infh, outbamfh, genome_rname_tid_map, tid_tx_map

def transcriptome_to_genome(genome_index,
                            transcript_feature_file, 
                            input_file, 
                            output_bam_file,
                            library_type,
                            input_is_sam):
    # setup and open files
    infh, outbamfh, genome_rname_tid_map, tid_tx_map = \
        _setup_and_open_files(genome_index, transcript_feature_file,
                              input_file, output_bam_file, library_type,
                              input_is_sam)
    # now convert BAM reads
    logging.debug("Converting transcriptome to genome BAM")
    for pe_reads in parse_pe_reads(infh):
        pairs, unpaired_reads = group_read_pairs(pe_reads)
        if len(pairs) > 0:
            # convert pairs
            pairs_dict = collections.OrderedDict()
            for r1,r2 in pairs:
                newr1 = convert_read(r1, tid_tx_map, genome_rname_tid_map)
                newr2 = convert_read(r2, tid_tx_map, genome_rname_tid_map)
                pair_reads(newr1, newr2)
                # key to identify independent alignments
                k = (newr1.rname, newr1.pos, newr1.aend, 
                     newr2.rname, newr2.pos, newr2.aend)
                if k not in pairs_dict:
                    pairs_dict[k] = (newr1, newr2)
            # compute number of alignment hits
            num_hits = len(pairs_dict)
            # write reads to BAM file
            for r1,r2 in pairs_dict.itervalues():
                # annotate multihits
                r1.tags = r1.tags + [("NH", num_hits)]
                r2.tags = r2.tags + [("NH", num_hits)]
                outbamfh.write(r1)
                outbamfh.write(r2)
        else:
            # convert unpaired reads
            unpaired_reads_dict = (collections.OrderedDict(),
                                   collections.OrderedDict())
            # compute number of alignment hits
            num_hits = [0,0]
            for rnum,reads in enumerate(unpaired_reads):
                for r in reads:
                    newr = convert_read(r, tid_tx_map, genome_rname_tid_map)
                    # key to identify independent alignments
                    k = (newr.rname, newr.pos, newr.aend)
                    if k not in unpaired_reads_dict[rnum]:
                        unpaired_reads_dict[rnum][k] = r
                        if not r.is_unmapped:
                            num_hits[rnum] += 1
            # compute number of alignment hits
            for rnum,reads_dict in enumerate(unpaired_reads_dict):                
                for r in reads_dict.itervalues():
                    # annotate multihits
                    r.tags = r.tags + [("NH", num_hits[rnum])]
                    outbamfh.write(r)
    outbamfh.close()
    infh.close()
    return config.JOB_SUCCESS

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("--library-type", dest="library_type",
                        default=LibraryTypes.FR_UNSTRANDED)
    parser.add_argument("--sam", dest="sam", action="store_true", 
                        default=False)
    parser.add_argument("genome_index")
    parser.add_argument("transcript_feature_file")
    parser.add_argument("input_sam_file")
    parser.add_argument("output_bam_file") 
    args = parser.parse_args()
    return transcriptome_to_genome(args.genome_index,
                                   args.transcript_feature_file, 
                                   args.input_sam_file, 
                                   args.output_bam_file,
                                   args.library_type,
                                   args.sam)

if __name__ == '__main__':
    sys.exit(main())
