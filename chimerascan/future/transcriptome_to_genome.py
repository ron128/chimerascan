'''
Created on Jun 6, 2012

@author: mkiyer
'''
import sys
import os
import logging
import argparse
import collections

from chimerascan import pysam
from chimerascan.lib import config
from chimerascan.lib.feature import TranscriptFeature
from chimerascan.lib.sam import copy_read

#def get_sam_header_from_bowtie2_index(index):
#    # extract sequence names and lengths from bowtie reference
#    args = [config.BOWTIE2_INSPECT_BIN, "-s", index]
#    p = subprocess.Popen(args, stdout=subprocess.PIPE)
#    stdoutdata = p.communicate()[0]
#    retcode = p.wait()
#    if retcode != 0:
#        return config.JOB_ERROR, {}
#    # parse index information and build SAM header
#    headerdict = {}
#    sqlist = []
#    for line in stdoutdata.split('\n'):
#        if not line:
#            continue
#        line = line.strip()
#        if not line:
#            continue
#        if not line.startswith("Sequence"):
#            continue
#        fields = line.strip().split('\t')
#        seqname = fields[1]
#        seqlen = int(fields[2])
#        sqlist.append({'LN': seqlen, 'SN': seqname})
#    headerdict['SQ'] = sqlist
#    return headerdict
#    # build SAM header from genome index
#    #if not check_executable(config.BOWTIE2_INSPECT_BIN):
#    #    logging.error("Cannot find bowtie2-inspect binary")
#    #    return config.JOB_ERROR
#    #bt2_genome_index = os.path.join(index_dir, config.GENOME_INDEX)
#    #headerdict = get_sam_header_from_bowtie2_index(bt2_genome_index)
#    # 

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
    for cigarcode, cigarbp in cigar:
        if ((cigarcode == CIGAR_M) or 
            (cigarcode == CIGAR_D) or
            (cigarcode == CIGAR_N)):
            # process the aligned cigar bp
            while cigarbp > (exon_size - toffset):
                # subtract remainder of exon from cigar bp
                cigarbp -= (exon_size - toffset)
                # add cigar for remainder of this exon
                newcigar.append((cigarcode, exon_size - toffset))
                # insert skip CIGAR operation for exon junction
                prev_estart, prev_eend = exons[eindex]
                eindex += 1
                if negstrand:
                    intron_size = (prev_estart - exons[eindex][1])
                else:
                    intron_size = (exons[eindex][0] - prev_eend)
                newcigar.append((CIGAR_N, intron_size))
                # advance to next exon
                testart += exon_size
                toffset = 0
                exon_size = exons[eindex][1] - exons[eindex][0]     
            # update offset into current exon
            toffset += cigarbp
        # can finish remaining cigarbp within current exon
        newcigar.append((cigarcode, cigarbp))
    if negstrand:
        # flip CIGAR list
        newcigar.reverse()        
    return newcigar

def transcriptome_to_genome(template_genome_bam_file,
                            transcript_feature_file, 
                            input_bam_file, 
                            output_bam_file):
    # copy genome BAM file header
    genomefh = pysam.Samfile(template_genome_bam_file, "rb")
    outbamfh = pysam.Samfile(output_bam_file, "wb", template=genomefh)
    genomefh.close()
    # convert transcriptome alignments to genomic positions    
    inbamfh = pysam.Samfile(input_bam_file, "rb")
    genome_rname_tid_map = dict((rname,i) for i,rname in enumerate(outbamfh.references))
    tx_rname_tid_map = dict((rname,i) for i,rname in enumerate(inbamfh.references))
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
    # now convert BAM reads
    logging.debug("Converting transcriptome to genome BAM")
    for r in inbamfh:
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
            newcigar = convert_cigar(r.cigar, negstrand, exons, eindex, 
                                     testart, toffset)            
            newr.rname = genome_tid
            newr.pos = newpos
            newr.cigar = newcigar
            # find genomic position of alignment
            if negstrand:                
                # flip is_reverse flag
                newr.is_reverse = (not newr.is_reverse)
        # remove mate information
        newr.is_proper_pair = False
        newr.mrnm = -1
        newr.mpos = 0
        newr.mate_is_unmapped = True
        newr.mate_is_reverse = True
        newr.isize = 0
#        if not r.mate_is_unmapped:
#            chrom, negstrand, exons = tid_tx_map[r.mrnm]
#            # TODO: remove assert statement
#            assert chrom in genome_rname_tid_map
#            # convert mate reference name and position 
#            mate_genome_tid = genome_rname_tid_map[chrom]
#            mate_pos = convert_pos(r.mpos, negstrand, exons)[0]
#            newr.mrnm = mate_genome_tid
#            newr.mpos = mate_pos
#            if negstrand:
#                # flip is_reverse flag
#                newr.mate_is_reverse = (not newr.mate_is_reverse)
#            else:
#                pass
        outbamfh.write(newr)            
    outbamfh.close()
    inbamfh.close()
    return config.JOB_SUCCESS

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("template_genome_bam_file")
    parser.add_argument("transcript_feature_file")
    parser.add_argument("input_bam_file")
    parser.add_argument("output_bam_file") 
    args = parser.parse_args()
    return transcriptome_to_genome(args.template_genome_bam_file,
                                   args.transcript_feature_file, 
                                   args.input_bam_file, 
                                   args.output_bam_file)

if __name__ == '__main__':
    sys.exit(main())
