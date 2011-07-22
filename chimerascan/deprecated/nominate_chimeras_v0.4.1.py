'''
Created on Jul 21, 2011

@author: mkiyer
'''

'''
Created on Jun 4, 2011

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
import os
import collections
import itertools
import operator

from chimerascan import pysam

from chimerascan.lib import config
from chimerascan.lib.chimera import Chimera, \
    DiscordantTags, DISCORDANT_TAG_NAME, \
    OrientationTags, ORIENTATION_TAG_NAME, \
    DiscordantRead, ChimeraTypes, ChimeraPartner
from chimerascan.lib.gene_to_genome import build_tid_tx_maps

def parse_pairs(bamfh):
    bam_iter = iter(bamfh)
    try:
        while True:
            r1 = bam_iter.next()
            r2 = bam_iter.next()
            yield r1,r2
    except StopIteration:
        pass

def parse_gene_chimeric_reads(bamfh):
    # create a dictionary structure to hold read pairs
    chimera_dict = collections.defaultdict(lambda: [])   
    for r1,r2 in parse_pairs(bamfh):
        #
        # TODO:
        # for now we are only going to deal with gene-gene
        # chimeras and leave other chimeras for study at a 
        # later time
        #
        dr1 = r1.opt(DISCORDANT_TAG_NAME)
        dr2 = r2.opt(DISCORDANT_TAG_NAME)
        if (dr1 != DiscordantTags.DISCORDANT_GENE or
            dr2 != DiscordantTags.DISCORDANT_GENE):            
            continue
        # organize key in 5' to 3' order
        or1 = r1.opt(ORIENTATION_TAG_NAME)
        or2 = r2.opt(ORIENTATION_TAG_NAME)
        assert or1 != or2
        if or1 == OrientationTags.FIVEPRIME:
            pair = (r1,r2)
        else:
            pair = (r2,r1)
        # store pertinent information in lightweight structure
        # convert to DiscordantRead objects
        r5p = DiscordantRead.from_read(pair[0])
        r3p = DiscordantRead.from_read(pair[1])
        # keep list of discordant chimeric reads
        chimera_dict[(r5p.tid, r3p.tid)].append((r5p,r3p))
    for key,pairs in chimera_dict.iteritems():
        rname1,rname2 = key
        yield rname1, rname2, pairs

def get_chimera_type(fiveprime_gene, threeprime_gene, gene_trees):
    """
    return tuple containing ChimeraType and distance 
    between 5' and 3' genes 
    """
    # get gene information
    chrom5p, start5p, end5p, strand1 = fiveprime_gene.chrom, fiveprime_gene.tx_start, fiveprime_gene.tx_end, fiveprime_gene.strand
    chrom3p, start3p, end3p, strand2 = threeprime_gene.chrom, threeprime_gene.tx_start, threeprime_gene.tx_end, threeprime_gene.strand
    # interchromosomal
    if chrom5p != chrom3p:
        return ChimeraTypes.INTERCHROMOSOMAL, None
    # orientation
    same_strand = strand1 == strand2
    # genes on same chromosome so check overlap
    is_overlapping = (start5p < end3p) and (start3p < end5p)            
    if is_overlapping:
        if not same_strand:
            if ((start5p <= start3p and strand1 == "+") or
                (start5p > start3p and strand1 == "-")):                    
                return (ChimeraTypes.OVERLAP_CONVERGE, 0)
            else:
                return (ChimeraTypes.OVERLAP_DIVERGE, 0)
        else:
            if ((start5p <= start3p and strand1 == "+") or
                (end5p >= end3p and strand1 == "-")):
                return (ChimeraTypes.OVERLAP_SAME, 0)
            else:
                return (ChimeraTypes.OVERLAP_COMPLEX, 0)
    # if code gets here then the genes are on the same chromosome but do not
    # overlap.  first calculate distance (minimum distance between genes)
    if start5p <= start3p:
        distance = start3p - end5p
        between_start,between_end = end5p,start3p
    else:
        distance = end3p - start5p
        between_start,between_end = end3p,start5p
    # check whether there are genes intervening between the
    # chimera candidates
    genes_between = []
    genes_between_same_strand = []
    for hit in gene_trees[chrom5p].find(between_start,
                                       between_end):
        if (hit.start > between_start and
            hit.end < between_end):             
            if hit.strand == strand1:
                genes_between_same_strand.append(hit)
            genes_between.append(hit)
            
    if same_strand:
        if len(genes_between_same_strand) == 0:
            return ChimeraTypes.READTHROUGH, distance
        else:
            return ChimeraTypes.INTRACHROMOSOMAL, distance
    else:
        # check for reads between neighboring genes    
        if len(genes_between) == 0:
            if ((start5p <= start3p and strand1 == "+") or
                (start5p > start3p and strand1 == "-")):                    
                return (ChimeraTypes.ADJ_CONVERGE, distance)
            elif ((start5p >= start3p and strand1 == "+") or
                  (start5p < start3p and strand1 == "-")):
                return (ChimeraTypes.ADJ_DIVERGE, distance)
            elif ((start5p <= start3p and strand1 == "+") or
                  (start5p > start3p and strand1 == "-")):
                return (ChimeraTypes.ADJ_SAME, distance)
            elif ((start5p >= start3p and strand1 == "+") or
                  (start5p < start3p and strand1 == '-')):
                return (ChimeraTypes.ADJ_COMPLEX, distance)
        else:
            return ChimeraTypes.INTRA_COMPLEX, distance    
    return ChimeraTypes.UNKNOWN, distance

def read_pairs_to_chimera(chimera_name, tid5p, tid3p, readpairs, 
                          tid_tx_map, genome_tx_trees, trim_bp):
    # get gene information
    tx5p = tid_tx_map[tid5p]
    tx3p = tid_tx_map[tid3p]
    # categorize chimera type
    chimera_type, distance = get_chimera_type(tx5p, tx3p, genome_tx_trees)
    # create chimera object
    c = Chimera()
    iter5p = itertools.imap(operator.itemgetter(0), readpairs)
    iter3p = itertools.imap(operator.itemgetter(1), readpairs)
    c.partner5p = ChimeraPartner.from_discordant_reads(iter5p, tx5p, trim_bp)
    c.partner3p = ChimeraPartner.from_discordant_reads(iter3p, tx3p, trim_bp)
    c.name = chimera_name
    c.chimera_type = chimera_type
    c.distance = distance
    # raw reads
    c.encomp_read_pairs = readpairs
    return c

def nominate_chimeras(index_dir, input_bam_file, output_file, trim_bp):
    logging.debug("Reading gene information")
    gene_file = os.path.join(index_dir, config.GENE_FEATURE_FILE)
    bamfh = pysam.Samfile(input_bam_file, "rb")
    # build a lookup table to get genomic intervals from transcripts
    tid_tx_map, genome_tx_trees = build_tid_tx_maps(bamfh, gene_file,
                                                    rname_prefix=config.GENE_REF_PREFIX)
    # group discordant read pairs by gene
    chimera_num = 0
    outfh = open(output_file, "w")
    logging.debug("Parsing discordant reads")
    for tid5p,tid3p,readpairs in parse_gene_chimeric_reads(bamfh):
        c = read_pairs_to_chimera("C%07d" % (chimera_num), tid5p, tid3p, 
                                  readpairs, tid_tx_map,
                                  genome_tx_trees, trim_bp)
        fields = c.to_list()
        chimera_num += 1
        print >>outfh, '\t'.join(map(str, fields))
    outfh.close()
    bamfh.close()

def main():
    from optparse import OptionParser
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = OptionParser("usage: %prog [options] <index> "
                          "<discordant_reads.srt.bedpe> <chimeras.txt>")
    parser.add_option("--trim", dest="trim", type="int", 
                      default=config.EXON_JUNCTION_TRIM_BP)
    options, args = parser.parse_args()
    index_dir = args[0]
    input_file = args[1]
    output_file = args[2]
    nominate_chimeras(index_dir, input_file, output_file, options.trim)


if __name__ == '__main__':
    main()
