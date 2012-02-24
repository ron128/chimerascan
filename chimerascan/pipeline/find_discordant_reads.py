'''
Created on Jun 2, 2011

@author: mkiyer
'''
import logging
import collections
import os

from chimerascan import pysam
from chimerascan.bx.cluster import ClusterTree

from chimerascan.lib import config
from chimerascan.lib.base import LibraryTypes
from chimerascan.lib.sam import parse_pe_reads, pair_reads, copy_read, select_best_mismatch_strata
from chimerascan.lib.gene_to_genome import build_rname_genome_map, \
    build_rname_cluster_map, gene_to_genome_pos
from chimerascan.lib.chimera import DiscordantTags, DISCORDANT_TAG_NAME, \
    OrientationTags, ORIENTATION_TAG_NAME

# globals
imin2 = lambda a,b: a if a <= b else b

def annotate_multihits(bamfh, reads, rname_genome_map, tid_rname_map):
    hits = set()
    any_unmapped = False
    for r in reads:
        if r.is_unmapped:
            any_unmapped = True
            continue
        # get reference name from tid
        rname = tid_rname_map[r.rname]
        # TODO: remove this after debugging
        assert rname in rname_genome_map
        # use the position that is smallest relative to genome
        left_rname, left_strand, left_pos = gene_to_genome_pos(rname, r.pos, rname_genome_map)
        right_rname, right_strand, right_pos = gene_to_genome_pos(rname, r.aend-1, rname_genome_map)
        pos = imin2(left_pos, right_pos)
        hits.add((left_rname, pos))
        #print r.qname, bamfh.getrname(r.rname), r.pos, bamfh.getrname(tid), pos  
    for i,r in enumerate(reads):
        # annotate reads with 'HI', and 'IH' tags
        r.tags = r.tags + [("HI",i), ("IH",len(reads)), ("NH", len(hits))]
    return any_unmapped

def map_reads_to_references(pe_reads, rname_cluster_map, tid_rname_map):
    """
    bin reads by transcript cluster and reference (tid)
    """
    refdict = collections.defaultdict(lambda: ([], []))
    genedict = collections.defaultdict(lambda: ([], []))
    for readnum, reads in enumerate(pe_reads):
        for r in reads:
            if r.is_unmapped:
                continue 
            # get reference name from tid
            rname = tid_rname_map[r.rname]
            # get cluster id
            if rname in rname_cluster_map:
                # add to cluster dict
                cluster_id = rname_cluster_map[rname]
                pairs = genedict[cluster_id]
                pairs[readnum].append(r)
            # add to reference dict
            pairs = refdict[r.rname]
            pairs[readnum].append(r)
    return refdict, genedict

def get_genome_orientation(r, library_type):
    if library_type == LibraryTypes.FR_FIRSTSTRAND:
        if r.is_read2:
            return OrientationTags.FIVEPRIME
        else:
            return OrientationTags.THREEPRIME
    elif library_type == LibraryTypes.FR_SECONDSTRAND:
        if r.is_read1:
            return OrientationTags.FIVEPRIME
        else:
            return OrientationTags.THREEPRIME
    return OrientationTags.NONE

def get_gene_orientation(r, library_type):
    if library_type == LibraryTypes.FR_UNSTRANDED:
        if r.is_reverse:
            return OrientationTags.THREEPRIME
        else:
            return OrientationTags.FIVEPRIME
    elif library_type == LibraryTypes.FR_FIRSTSTRAND:
        if r.is_read2:
            return OrientationTags.FIVEPRIME
        else:
            return OrientationTags.THREEPRIME
    elif library_type == LibraryTypes.FR_SECONDSTRAND:
        if r.is_read1:
            return OrientationTags.FIVEPRIME
        else:
            return OrientationTags.THREEPRIME
    logging.error("Unknown library type %s, aborting" % (library_type))
    assert False

def classify_unpaired_reads(reads, library_type):
    gene_hits_5p = []
    gene_hits_3p = []
    for r in reads:
        # this alignment is to a transcript (gene), so need
        # to determine whether it is 5' or 3'
        orientation = get_gene_orientation(r, library_type)
        if orientation == OrientationTags.FIVEPRIME:
            gene_hits_5p.append(r)
        else:
            gene_hits_3p.append(r)
        # add a tag to the sam file describing the read orientation and
        # that it is discordant
        r.tags = r.tags + [(DISCORDANT_TAG_NAME, DiscordantTags.DISCORDANT_GENE),
                           (ORIENTATION_TAG_NAME, orientation)]                               
    return gene_hits_5p, gene_hits_3p

def find_discordant_pairs(pe_reads, library_type):
    """
    iterate through combinations of read1/read2 to predict valid 
    discordant read pairs
    """
    # classify the reads as 5' or 3' gene alignments or genome alignments
    r1_5p_gene_hits, r1_3p_gene_hits = \
        classify_unpaired_reads(pe_reads[0], library_type)
    r2_5p_gene_hits, r2_3p_gene_hits = \
        classify_unpaired_reads(pe_reads[1], library_type)
    # pair 5' and 3' gene alignments
    gene_pairs = []
    combos = [(r1_5p_gene_hits,r2_3p_gene_hits),
              (r1_3p_gene_hits,r2_5p_gene_hits)]
    for r1_list,r2_list in combos:
        for r1 in r1_list:
            for r2 in r2_list:
                cr1 = copy_read(r1)
                cr2 = copy_read(r2)
                pair_reads(cr1,cr2)
                gene_pairs.append((cr1,cr2))
    return gene_pairs

def classify_read_pairs(pe_reads, max_isize,
                        library_type,
                        rname_cluster_map,
                        tid_rname_map):
    """
    examines all the alignments of a single fragment and tries to find ways
    to pair reads together.
    
    annotates all read pairs with an integer tag corresponding to a value
    in the DiscordantTags class
    
    returns a tuple with the following lists:
    1) pairs (r1,r2) aligning to genes (pairs may be discordant)
    """
    # to satisfy library type reads must either be on 
    # same strand or opposite strands
    concordant_tx_pairs = []
    discordant_tx_pairs = []
    concordant_gene_pairs = []
    discordant_gene_pairs = []
    # 
    # first, try to pair reads that map to the same transcript, or to the
    # genome within the insert size range
    #
    same_strand = LibraryTypes.same_strand(library_type)
    refdict,clusterdict = map_reads_to_references(pe_reads, rname_cluster_map, tid_rname_map)
    found_pair = False
    for tid, tid_pe_reads in refdict.iteritems():
        # check if there are alignments involving both reads in a pair
        if len(tid_pe_reads[0]) == 0 or len(tid_pe_reads[1]) == 0:
            # no paired alignments exist at this reference
            continue
        # check if there are alignments involving both reads in a pair        
        for r1 in tid_pe_reads[0]:
            for r2 in tid_pe_reads[1]:
                # these reads can be paired
                # TODO: insist on being within an expectedc insert size distance?
                found_pair = True
                cr1 = copy_read(r1)
                cr2 = copy_read(r2)                    
                # this is a hit to same transcript (gene)
                # read strands must agree with library type
                # reads are concordant if strand comparison is correct
                strand_match = (same_strand == (r1.is_reverse == r2.is_reverse))
                if strand_match:
                    tags = [(DISCORDANT_TAG_NAME, DiscordantTags.CONCORDANT_TX)]
                    concordant_tx_pairs.append((cr1,cr2))
                else:
                    # hit to same gene with wrong strand, which
                    # could happen in certain wacky cases
                    tags = [(DISCORDANT_TAG_NAME, DiscordantTags.DISCORDANT_STRAND_TX)]
                    discordant_tx_pairs.append((cr1,cr2))
                pair_reads(cr1,cr2,tags)
    # at this point, if we have not been able to find a suitable way
    # to pair the reads, then search within the transcript cluster
    if not found_pair:
        for cluster_pe_reads in clusterdict.itervalues():
            # check if there are alignments involving both reads in a pair
            if len(cluster_pe_reads[0]) == 0 or len(cluster_pe_reads[1]) == 0:
                # no paired alignments in this transcript cluster            
                continue
            for r1 in cluster_pe_reads[0]:
                for r2 in cluster_pe_reads[1]:
                    # check strand compatibility
                    strand_match = (same_strand == (r1.is_reverse == r2.is_reverse))
                    # these reads can be paired
                    found_pair = True
                    cr1 = copy_read(r1)
                    cr2 = copy_read(r2)                    
                    if strand_match:
                        tags = [(DISCORDANT_TAG_NAME, DiscordantTags.CONCORDANT_GENE)]
                        concordant_gene_pairs.append((cr1,cr2))
                    else:
                        tags = [(DISCORDANT_TAG_NAME, DiscordantTags.DISCORDANT_STRAND_GENE)]
                        discordant_gene_pairs.append((cr1,cr2))
                    pair_reads(cr1,cr2,tags)
    # at this point, we have tried all combinations.  if any paired reads
    # are concordant then return them without considering discordant reads 
    gene_pairs = []
    if len(concordant_tx_pairs) > 0:
        gene_pairs = concordant_tx_pairs
    elif len(concordant_gene_pairs) > 0:
        gene_pairs = concordant_gene_pairs
    if len(gene_pairs) > 0:
        return gene_pairs
    # if no concordant reads in transcripts or genome, return any
    # discordant reads that may violate strand requirements but still
    # remain colocalized on the same gene/chromosome
    gene_pairs = []
    if len(discordant_tx_pairs) > 0:
        gene_pairs = discordant_tx_pairs
    elif len(discordant_gene_pairs) > 0:
        gene_pairs = discordant_gene_pairs    
    if len(gene_pairs) > 0:
        return gene_pairs
    #
    # at this point, no read pairings were found so the read is 
    # assumed to be discordant.  
    #
    # TODO: now that we know that the reads are discordant, no reason
    # to keep all the mappings hanging around if there is a small subset
    # with a small number of mismatches.  is this the right thing to do
    # here?
    # 
    pe_reads = (select_best_mismatch_strata(pe_reads[0]),
                select_best_mismatch_strata(pe_reads[1]))
    #
    # now we can create all valid combinations of read1/read2 as putative 
    # discordant read pairs 
    #    
    gene_pairs = find_discordant_pairs(pe_reads, library_type)
    return gene_pairs

def write_pe_reads(bamfh, pe_reads):
    for reads in pe_reads:
        for r in reads:
            bamfh.write(r)

def write_pairs(bamfh, pairs):
    for r1,r2 in pairs:
        bamfh.write(r1)
        bamfh.write(r2)

def find_discordant_fragments(input_bam_file, gene_paired_bam_file,
                              unmapped_bam_file, index_dir, max_isize, 
                              library_type):
    """
    parses BAM file and categorizes reads into several groups:
    - concordant
    - discordant within gene (splicing isoforms)
    - discordant between different genes (chimeras)
    - discordant genome alignments (unannotated)
    """
    logging.info("Finding discordant read pair combinations")
    logging.debug("\tInput file: %s" % (input_bam_file))
    logging.debug("\tMax insert size: '%d'" % (max_isize))
    logging.debug("\tLibrary type: '%s'" % (library_type))
    logging.debug("\tGene paired file: %s" % (gene_paired_bam_file))
    logging.debug("\tUnmapped file: %s" % (unmapped_bam_file))
    # setup input and output files
    bamfh = pysam.Samfile(input_bam_file, "rb")
    genefh = pysam.Samfile(gene_paired_bam_file, "wb", template=bamfh)
    unmappedfh = pysam.Samfile(unmapped_bam_file, "wb", template=bamfh)
    gene_file = os.path.join(index_dir, config.GENE_FEATURE_FILE)
    # build a lookup table to get all the overlapping transcripts given a
    # transcript 'tid'
    rname_cluster_map = build_rname_cluster_map(open(gene_file), 
                                                rname_prefix=config.GENE_REF_PREFIX)
    # build a lookup table to get genome coordinates from transcript 
    # coordinates
    rname_genome_map = build_rname_genome_map(open(gene_file), 
                                              rname_prefix=config.GENE_REF_PREFIX)
    # get tid -> rname mappings
    tid_rname_map = bamfh.references
    for pe_reads in parse_pe_reads(bamfh):
        # add hit index and number of multimaps information to read tags
        # this function also checks for unmapped reads
        any_unmapped = False
        for reads in pe_reads:
            any_unmapped = (any_unmapped or 
                            annotate_multihits(bamfh, reads, 
                                               rname_genome_map, 
                                               tid_rname_map))
        if any_unmapped:
            # write to output as discordant reads and continue to 
            # next fragment
            write_pe_reads(unmappedfh, pe_reads)
            continue
        # examine all read pairing combinations and rule out invalid 
        # pairings.  this returns gene pairs and genome pairs
        gene_pairs = classify_read_pairs(pe_reads, max_isize, library_type,
                                         rname_cluster_map, tid_rname_map)
        write_pairs(genefh, gene_pairs)
    genefh.close()
    unmappedfh.close()
    bamfh.close()  
    logging.info("Finished pairing reads")

def main():
    from optparse import OptionParser
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = OptionParser("usage: %prog [options] <index> <in.bam> "
                          "<gene_paired.bam> <unmapped.bam>")
    parser.add_option('--max-fragment-length', dest="max_fragment_length", 
                      type="int", default=1000)
    parser.add_option('--library', dest="library_type", 
                      default=LibraryTypes.FR_UNSTRANDED)
    options, args = parser.parse_args()    
    index_dir = args[0]
    input_bam_file = args[1]
    gene_paired_bam_file = args[2]
    unmapped_bam_file = args[3]
    find_discordant_fragments(input_bam_file, gene_paired_bam_file,
                              unmapped_bam_file, index_dir,
                              max_isize=options.max_fragment_length,
                              library_type=options.library_type)

if __name__ == '__main__':
    main()