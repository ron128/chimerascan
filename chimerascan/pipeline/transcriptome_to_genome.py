'''
Created on Jun 6, 2012

@author: mkiyer
'''
import sys
import logging
import argparse
import collections
import subprocess

import pysam

import chimerascan
from chimerascan.lib import config
from chimerascan.lib.seq import DNA_reverse_complement
from chimerascan.lib.base import check_executable, LibraryTypes
from chimerascan.lib.feature import TranscriptFeature
from chimerascan.lib.sam import copy_read, parse_pe_reads, \
    group_read_pairs, pair_reads, REF_ADVANCING_CIGAR_CODES, CIGAR_N

def get_references_from_bowtie2_index(index):
    # extract sequence names and lengths from bowtie reference
    args = [config.BOWTIE2_INSPECT_BIN, "-s", index]
    p = subprocess.Popen(args, stdout=subprocess.PIPE)
    stdoutdata = p.communicate()[0]
    retcode = p.wait()
    if retcode != 0:
        return config.JOB_ERROR, {}
    # parse index information and format as SAM header
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
        sqlist.append((seqname, seqlen))
    return sqlist

def reverse_complement_MD_tag(md):
    digits = []
    mdops = []
    for i in xrange(len(md)):
        if md[i].isdigit():
            digits.append(md[i])
        else:
            if len(digits) > 0:
                mdops.append(''.join(digits))
                digits = []
            if md[i] == '^':
                mdops.append(md[i])
            elif md[i].isalpha():
                mdops.append(DNA_reverse_complement(md[i]))
    if len(digits) > 0:
        mdops.append(''.join(digits))
    return ''.join(mdops[::-1])


_library_type_strand_map = ({LibraryTypes.FR_UNSTRANDED: {True: '.', False: '.'},
                             LibraryTypes.FR_FIRSTSTRAND: {True: '+', False: '-'},
                             LibraryTypes.FR_SECONDSTRAND: {True: '-', False: '+'}},
                            {LibraryTypes.FR_UNSTRANDED: {True: '.', False: '.'},
                             LibraryTypes.FR_FIRSTSTRAND: {True: '-', False: '+'},
                             LibraryTypes.FR_SECONDSTRAND: {True: '+', False: '-'}})
def get_read_strand(is_read2, is_reverse, negstrand, library_type):
    if library_type == LibraryTypes.FR_UNSTRANDED:
        if negstrand:
            strand = '-'
        else:
            strand = '+'
    else:
        rnum = int(is_read2)
        strand = _library_type_strand_map[rnum][library_type][is_reverse]
    return strand 

def convert_pos(pos, negstrand, exons):
    """
    find position along transcript
    """
    eindex = 0
    testart = 0     
    exon_size = exons[eindex][1] - exons[eindex][0]
    while pos >= (testart + exon_size):
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
    alen = 0
    spliced = False
    for cigarcode, cigarbp in cigar:
        if cigarcode in REF_ADVANCING_CIGAR_CODES:
            # process the aligned cigar bp
            while cigarbp > (exon_size - toffset):
                # subtract remainder of exon from cigar bp
                cigarbp -= (exon_size - toffset)
                # add cigar for remainder of this exon
                newcigar.append((cigarcode, exon_size - toffset))
                alen += (exon_size - toffset)                
                # insert skip CIGAR operation for exon junction
                prev_estart, prev_eend = exons[eindex]
                eindex += 1
                if negstrand:
                    intron_size = (prev_estart - exons[eindex][1])
                else:
                    intron_size = (exons[eindex][0] - prev_eend)
                newcigar.append((CIGAR_N, intron_size))
                alen += intron_size
                spliced = True     
                # advance to next exon
                testart += exon_size
                toffset = 0
                exon_size = exons[eindex][1] - exons[eindex][0]     
            # update offset into current exon
            toffset += cigarbp
            alen += cigarbp
        # can finish remaining cigarbp within current exon
        newcigar.append((cigarcode, cigarbp))
    if negstrand:
        # flip CIGAR list
        newcigar.reverse() 
    return newcigar, alen, spliced

def convert_read(r, transcript_tid_map, library_type):
    if r.is_unmapped:
        # return copy of original read
        return copy_read(r)
    # copy and modify tags
    tagdict = collections.OrderedDict(r.tags)
    if 'XS' in tagdict:
        del tagdict['XS']
    if 'NH' in tagdict:
        del tagdict['NH']
    # convert transcript reference to genome
    genome_tid, negstrand, exons = transcript_tid_map[r.tid]
    # find genomic start position of transcript
    newpos, eindex, testart, toffset = convert_pos(r.pos, negstrand, exons)
    # parse and convert transcript cigar string
    newcigar, alen, spliced = \
        convert_cigar(r.cigar, negstrand, exons, 
                      eindex, testart, toffset)            
    if negstrand:
        # set position to left end of transcript
        newpos = newpos - alen + 1            
        # flip is_reverse flag
        is_reverse = (not r.is_reverse)
        # reverse complement seq and quals
        seq = DNA_reverse_complement(r.seq)
        qual = None if r.qual is None else r.qual[::-1]
        # flip MD tag
        if 'MD' in tagdict:
            tagdict['MD'] = reverse_complement_MD_tag(tagdict['MD'])
    else:
        is_reverse = r.is_reverse
        seq = r.seq
        qual = r.qual
    # add XS tag
    strand = get_read_strand(r.is_read2, is_reverse, negstrand, library_type)
    tagdict['XS'] = strand
    # create copy of read
    a = pysam.AlignedRead()
    a.qname = r.qname
    a.flag = r.flag
    a.seq = seq
    a.qual = qual
    a.is_reverse = is_reverse
    a.tid = genome_tid
    a.pos = newpos
    a.cigar = newcigar
    a.mapq = r.mapq
    a.rnext = r.rnext
    a.pnext = r.pnext
    a.tlen = r.tlen
    a.tags = tuple(tagdict.iteritems())
    return a

def convert_read_pairs(pairs, transcript_tid_map, library_type):
    # convert pairs
    pairs_dict = collections.OrderedDict()
    for r1,r2 in pairs:
        newr1 = convert_read(r1, transcript_tid_map, library_type)
        newr2 = convert_read(r2, transcript_tid_map, library_type)
        pair_reads(newr1, newr2)
        # key to identify independent alignments
        k = (newr1.tid, newr1.pos, newr1.aend, 
             newr2.tid, newr2.pos, newr2.aend)
        if k not in pairs_dict:
            pairs_dict[k] = (newr1, newr2)
    # compute number of alignment hits
    num_hits = len(pairs_dict)
    # write reads to BAM file
    for r1,r2 in pairs_dict.itervalues():
        tagdict1 = collections.OrderedDict(r1.tags)
        tagdict2 = collections.OrderedDict(r2.tags)
        # annotate multihits
        tagdict1['NH'] = num_hits
        tagdict2['NH'] = num_hits
        # write
        r1.tags = tagdict1.items() 
        r2.tags = tagdict2.items()
        yield r1,r2

def convert_unpaired_reads(pe_reads, transcript_tid_map, library_type):
    # convert unpaired reads
    unpaired_reads_dict = (collections.OrderedDict(),
                           collections.OrderedDict())
    # compute number of alignment hits
    mate_num_hits = [0, 0]
    for rnum,reads in enumerate(pe_reads):
        for r in reads:
            newr = convert_read(r, transcript_tid_map, library_type)
            # key to identify independent alignments
            k = (newr.tid, newr.pos, newr.aend)
            if k not in unpaired_reads_dict[rnum]:
                unpaired_reads_dict[rnum][k] = newr
                if not r.is_unmapped:
                    mate_num_hits[rnum] += 1
    # compute number of alignment hits
    for rnum,reads_dict in enumerate(unpaired_reads_dict):                
        for r in reads_dict.itervalues():
            tagdict = collections.OrderedDict(r.tags)
            # annotate multihits
            tagdict['NH'] = mate_num_hits[rnum]
            # annotate multihits
            r.tags = tuple(tagdict.iteritems())
            yield r

def _setup_and_open_files(genome_index, transcripts,
                          input_file, output_file, 
                          library_type, input_sam, 
                          output_sam):
    # create SAM header from genome index
    logging.debug("Creating genome SAM header")
    if not check_executable(config.BOWTIE2_INSPECT_BIN):
        logging.error("Cannot find bowtie2-inspect binary")
        return config.JOB_ERROR
    # get references/lengths from bowtie2
    ref_list = get_references_from_bowtie2_index(genome_index)
    # open input BAM file and add to header
    if input_sam:
        mode = "r"
    else:
        mode = "rb"
    infh = pysam.Samfile(input_file, mode)
    header_dict = dict(infh.header)
    header_dict['SQ'] = [{'SN': seqname, 'LN': seqlen} for seqname,seqlen in ref_list]
    # open output BAM file with new header
    if output_sam:
        mode = "wh"
    else:
        mode = "wb"
    outfh = pysam.Samfile(output_file, mode, header=header_dict)
    # setup reference name mappings
    genome_rname_tid_map = dict((rname,i) for i,rname in enumerate(outfh.references))    
    transcriptome_rname_tid_map = dict((rname,i) for i,rname in enumerate(infh.references))
    # read transcript feature and prepare data structure for conversion
    logging.debug("Creating transcript to genome map")
    transcript_tid_map = {}
    for t in transcripts:
        exons = [(start, end) for start, end in t.exons]
        negstrand = True if t.strand == "-" else False
        if negstrand:
            exons.reverse()
        transcript_tid = transcriptome_rname_tid_map[str(t.tx_id)]
        genome_tid = genome_rname_tid_map[t.chrom]
        transcript_tid_map[transcript_tid] = (genome_tid, negstrand, exons)        
    return infh, outfh, transcript_tid_map

def transcriptome_to_genome(genome_index,
                            transcripts, 
                            input_file, 
                            output_file,
                            library_type,
                            input_sam,
                            output_sam):
    # setup and open files
    infh, outfh, transcript_tid_map = \
        _setup_and_open_files(genome_index, transcripts,
                              input_file, output_file, library_type,
                              input_sam, output_sam)
    # now convert BAM reads
    logging.debug("Converting transcriptome to genome BAM")
    num_paired_frags = 0
    num_unpaired_frags = 0
    for pe_reads in parse_pe_reads(infh):
        pairs, unpaired_reads = group_read_pairs(pe_reads)
        if len(pairs) > 0:
            num_paired_frags += 1
            # convert pairs
            for r1,r2 in convert_read_pairs(pairs, transcript_tid_map, 
                                            library_type):
                outfh.write(r1)
                outfh.write(r2)
        else:
            num_unpaired_frags += 1
            for r in convert_unpaired_reads(unpaired_reads, 
                                            transcript_tid_map, 
                                            library_type):
                outfh.write(r)
    logging.debug("Paired fragments: %d" % (num_paired_frags))
    logging.debug("Unpaired fragments: %d" % (num_unpaired_frags))
    outfh.close()
    infh.close()
    return config.JOB_SUCCESS

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("--library-type", dest="library_type",
                        default=LibraryTypes.FR_UNSTRANDED)
    parser.add_argument("--input-sam", dest="input_sam", action="store_true", 
                        default=False)
    parser.add_argument("--output-sam", dest="output_sam", action="store_true", 
                        default=False)
    parser.add_argument("genome_index")
    parser.add_argument("transcript_feature_file")
    parser.add_argument("input_sam_file")
    parser.add_argument("output_sam_file") 
    args = parser.parse_args()
    # read transcript features
    logging.debug("Reading transcript features")
    transcripts = list(TranscriptFeature.parse(open(args.transcript_feature_file)))
    return transcriptome_to_genome(args.genome_index,
                                   transcripts,
                                   args.input_sam_file, 
                                   args.output_sam_file,
                                   args.library_type,
                                   args.input_sam,
                                   args.output_sam)

if __name__ == '__main__':
    sys.exit(main())
