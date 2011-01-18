'''
Created on Jan 9, 2011

@author: mkiyer
'''
import collections
import logging
import os

import pysam

import config
from base import parse_library_type

def parse_pe_sam_file(bamfh):
    pe_reads = ([], [])
    # reads must be binned by qname, mate, hit, and segment
    # so initialize to mate 0, hit 0, segment 0
    num_reads = 0
    prev_qname = None
    for read in bamfh:
        # get read attributes
        qname = read.qname
        mate = 0 if read.is_read1 else 1
        # get hit/segment/mapping tags
        num_split_partitions = read.opt('NH')
        partition_ind = read.opt('XH')
        num_splits = read.opt('XN')
        split_ind = read.opt('XI')
        num_mappings = read.opt('IH')
        mapping_ind = read.opt('HI')        
        # if query name changes we have completely finished
        # the fragment and can reset the read data
        if num_reads > 0 and qname != prev_qname:
            yield pe_reads
            # reset state variables
            pe_reads = ([], [])
            num_reads = 0
        prev_qname = qname
        # initialize mate hits
        if len(pe_reads[mate]) == 0:
            pe_reads[mate].extend([list() for x in xrange(num_split_partitions)])
        mate_reads = pe_reads[mate]
        # initialize hit segments
        if len(mate_reads[partition_ind]) == 0:
            mate_reads[partition_ind].extend([list() for x in xrange(num_splits)])
        split_reads = mate_reads[partition_ind][split_ind]
        # initialize segment mappings
        if len(split_reads) == 0:
            split_reads.extend([None for x in xrange(num_mappings)])
        # add segment to hit/mate/read
        split_reads[mapping_ind] = read
        num_reads += 1
    if num_reads > 0:
        yield pe_reads

def map_reads_to_references(pe_reads):
    # bin reads by reference name to find reads that pairs
    # to the same gene/chromosome
    ref_dict = collections.defaultdict(lambda: ([], []))
    for mate, mate_hits in enumerate(pe_reads):
        # matching paired-end reads cannot have splits
        # to multiple references, and if multiple split
        # partitions were found it suggests split read
        # mapping occurred
        num_split_partitions = len(mate_hits)
        if num_split_partitions > 1:
            continue
        # reads with >1 split cannot be paired successfully
        # so do not add to reference dict
        if len(mate_hits[0]) > 1:
            continue            
        # this read has a single partition of splits and is
        # not split into multiple reads
        split_reads = mate_hits[0][0]
        for r in split_reads:
            if r.is_unmapped:
                continue 
            # add to reference dict
            mate_pairs = ref_dict[r.rname]
            mate_pairs[mate].append(r)
    return ref_dict

def find_concordant_pairs(ref_dict, min_isize, max_isize,
                          library_type):    
    same_strand = (library_type[0] == library_type[1])
    # check for mapping to same gene within insert size range
    concordant_pairs = []
    for rname, mate_pairs in ref_dict.iteritems():
        # both pairs must map to same reference
        if len(mate_pairs[0]) == 0 or len(mate_pairs[1]) == 0:
            continue
        # ensure distance is within insert size range
        # and strandedness matches library type
        for r1 in mate_pairs[0]:
            for r2 in mate_pairs[1]:
                # check insert size                                         
                if r1.pos > r2.pos:
                    isize = r1.aend - r2.pos
                else:
                    isize = r2.aend - r1.pos
                if isize < min_isize or isize > max_isize:
                    continue                
                # read strands must agree with library type
                if same_strand != (r1.is_reverse == r2.is_reverse):
                    continue                        
                # this is a concordant read pair
                concordant_pairs.append((r1, r2))
    return concordant_pairs

def get_insert_size():
    pass

def select_best_pairs(mate1_reads, mate2_reads,
                      min_fragment_length,
                      max_fragment_length,
                      library_type):
    pass

def pair_reads(r1, r2, add_tags=None, keep_tags=None):
    '''
    fill in paired-end fields in SAM record
    '''
    if keep_tags is None:
        keep_tags = []
    if add_tags is None:
        add_tags = []
    # convert read1 to paired-end
    r1.is_paired = True
    r1.is_proper_pair = True
    r1.is_read1 = True
    r1.mate_is_reverse = r2.is_reverse
    r1.mate_is_unmapped = r2.is_unmapped
    r1.mpos = r2.pos
    r1.mrnm = r2.rname
    # convert read2 to paired-end        
    r2.is_paired = True
    r2.is_proper_pair = True
    r2.is_read2 = True
    r2.mate_is_reverse = r1.is_reverse
    r2.mate_is_unmapped = r1.is_unmapped
    r2.mpos = r1.pos
    r2.mrnm = r1.rname
    # compute insert size
    if r1.pos > r2.pos:
        isize = r1.aend - r2.pos
    else:
        isize = r2.aend - r1.pos
    r1.isize = isize
    r2.isize = isize
    # update tags
    r1_tags = []
    r2_tags = []
    for tagname in keep_tags:
        r1_tags.append((tagname, r1.opt(tagname)))
        r2_tags.append((tagname, r2.opt(tagname)))        
    r1_tags.extend(add_tags)
    r2_tags.extend(add_tags)    
    r1.tags = r1_tags
    r2.tags = r2_tags
            
def merge_read_pairs(bamfh, output_bamfh, min_isize, max_isize, library_type):    
    # setup debugging logging messages
    debug_count = 0
    debug_every = 1e5
    debug_next = debug_every    
    num_paired = 0
    num_unpaired = 0
    num_fragments = 0
    for pe_reads in parse_pe_sam_file(bamfh):        
        ref_dict = map_reads_to_references(pe_reads)
        concordant_pairs = find_concordant_pairs(ref_dict, min_isize, max_isize,
                                                 library_type)        
        if len(concordant_pairs) > 0:
            for r1,r2 in concordant_pairs:
                pair_reads(r1, r2, 
                           keep_tags=('NM', 'MD'), 
                           add_tags=(('NH', len(concordant_pairs)),))
                output_bamfh.write(r1)
                output_bamfh.write(r2)
            # TODO: filter to select best pairs (fewest mismatches, insert size, etc)
            num_paired += 1
        else:
            # write unpaired reads to unpaired BAM file
            for mate_hits in pe_reads:
                for partitions in mate_hits:
                    for split_reads in partitions:
                        for r in split_reads:
                            output_bamfh.write(r)
            num_unpaired += 1
        num_fragments += 1
        # progress log
        debug_count += 1
        if debug_count == debug_next:
            debug_next += debug_every
            logging.info("Total read pairs: %d" % (num_fragments))
            logging.info("Paired reads: %d" % (num_paired))
            logging.info("Unpaired_reads: %d" % (num_unpaired))
    logging.info("Total read pairs: %d" % (num_fragments))
    logging.info("Paired reads: %d" % (num_paired))
    logging.info("Unpaired_reads: %d" % (num_unpaired))

def main():
    from optparse import OptionParser
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = OptionParser("usage: %prog [options] <in.bam> <out.bam>")
    parser.add_option('--min-fragment-length', dest="min_fragment_length", 
                      type="int", default=50)
    parser.add_option('--max-fragment-length', dest="max_fragment_length", 
                      type="int", default=1000)
    parser.add_option('--library', dest="library_type", default="fr")
    #parser.add_option('--unpaired-bam', dest="unpaired_bam_file", default=None)    
    options, args = parser.parse_args()
    input_bam_file = args[0]
    output_bam_file = args[1]
    logging.info("Merging read pairs")
    logging.debug("Input file: %s" % (input_bam_file))
    logging.debug("Output file: %s" % (output_bam_file))
    logging.debug("Library type: '%s'" % (options.library_type))
    library_type = parse_library_type(options.library_type)
    bamfh = pysam.Samfile(input_bam_file, "rb")
    outfh = pysam.Samfile(output_bam_file, "wb", template=bamfh)
    #outfh = pysam.Samfile("-", "w", template=bamfh)
    merge_read_pairs(bamfh, outfh, 
                     options.min_fragment_length,
                     options.max_fragment_length,
                     library_type)
    logging.info("Paired-end merging completed")
    
if __name__ == '__main__':
    main()
