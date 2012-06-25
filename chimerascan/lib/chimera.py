'''
Created on Jun 3, 2011

@author: mkiyer
'''
import logging
import collections

from base import LibraryTypes

STRAND_TAG = "XS"
STRAND_POS = "+"
STRAND_NEG = "-"
STRAND_NONE = "."

DISCORDANT_TAG_NAME = "XC"
class DiscordantTags(object):
    CONCORDANT_TX = 0
    DISCORDANT_STRAND_TX = 1
    CONCORDANT_GENE = 2
    DISCORDANT_STRAND_GENE = 3
    CONCORDANT_GENOME = 4
    DISCORDANT_STRAND_GENOME = 5
    DISCORDANT_GENE = 9
    DISCORDANT_GENOME = 17

ORIENTATION_TAG = "XD"
ORIENTATION_NONE = 0
ORIENTATION_5P = 1
ORIENTATION_3P = 2

def cmp_orientation(a,b):
    if (a == ORIENTATION_NONE) or (b == ORIENTATION_NONE):
        return True
    return (a != b)

def get_orientation(r, library_type):
    if library_type == LibraryTypes.FR_UNSTRANDED:
        if r.is_reverse:
            return ORIENTATION_3P
        else:
            return ORIENTATION_5P
    elif library_type == LibraryTypes.FR_FIRSTSTRAND:
        if r.is_read2:
            return ORIENTATION_5P
        else:
            return ORIENTATION_3P
    elif library_type == LibraryTypes.FR_SECONDSTRAND:
        if r.is_read1:
            return ORIENTATION_5P
        else:
            return ORIENTATION_3P
    logging.error("Unknown library type %s, aborting" % (library_type))
    assert False


DISCORDANT_CLUSTER_TAG = "XE"
DiscordantCluster = collections.namedtuple('DiscordantCluster', 
                                           ('rname', 'start', 'end', 
                                            'cluster_id', 'strand',
                                            'orientation', 'qnames',
                                            'concordant_frags'))

DiscordantClusterPair = collections.namedtuple('DiscordantClusterPair',
                                               ('pair_id', 'id5p', 'id3p', 'qnames'))

def parse_discordant_cluster_pair_file(line_iter):
    for line in line_iter:
        fields = line.strip().split('\t')
        pair_id = int(fields[0])
        id5p = int(fields[1])
        id3p = int(fields[2])
        qnames = fields[3].split(',')
        yield DiscordantClusterPair(pair_id, id5p, id3p, qnames)

class Chimera(object):
    _fields = ('rname5p', 'start5p', 'end5p',
               'rname3p', 'start3p', 'end3p',
               'chimera_id', 'num_frags',
               'strand5p', 'strand3p',
               'chimera_type', 'distance',
               'num_discordant_frags_5p',
               'num_discordant_frags_3p',
               'num_concordant_frags_5p',
               'num_concordant_frags_3p',
               'biotypes_5p', 'biotypes_3p',
               'genes_5p', 'genes_3p',
               'transcripts_5p', 'transcripts_3p')
    
    def __str__(self):
        fields = [self.rname5p, self.start5p, self.end5p,
                  self.rname3p, self.start3p, self.end3p,
                  self.chimera_id, self.num_frags,
                  self.strand5p, self.strand3p,
                  self.chimera_type, self.distance,
                  self.num_discordant_frags_5p,
                  self.num_discordant_frags_3p,
                  self.num_concordant_frags_5p,
                  self.num_concordant_frags_3p,
                  ','.join(self.biotypes_5p), 
                  ','.join(self.biotypes_3p),
                  ','.join(self.genes_5p),
                  ','.join(self.genes_3p),
                  ','.join(self.transcripts_5p),
                  ','.join(self.transcripts_3p)]
        return '\t'.join(map(str, fields))

    @staticmethod
    def from_string(line):
        fields = line.strip().split('\t')
        # format transcript information
        c = Chimera()
        c.rname5p = fields[0]
        c.start5p = int(fields[1])
        c.end5p = int(fields[2])
        c.rname3p = fields[3]
        c.start3p = int(fields[4])
        c.end3p = int(fields[5])
        c.chimera_id = fields[6]
        c.num_frags = int(fields[7])
        c.strand5p = fields[8]
        c.strand3p = fields[9]
        c.chimera_type = fields[10]
        c.distance = int(fields[11])
        c.num_discordant_frags_5p = int(fields[12])
        c.num_discordant_frags_3p = int(fields[13])
        c.num_concordant_frags_5p = int(fields[14])
        c.num_concordant_frags_3p = int(fields[15])
        c.biotypes_5p = fields[16].split(',')
        c.biotypes_3p = fields[17].split(',')
        c.genes_5p = fields[18].split(',')
        c.genes_3p = fields[19].split(',')
        c.transcripts_5p = fields[20].split(',')
        c.transcripts_3p = fields[21].split(',')
        return c
    
    @staticmethod
    def parse(line_iter):
        for line in line_iter:
            if line.startswith('#'):
                continue
            yield Chimera.from_string(line)

class ChimeraTypes(object):
    """
    chimera classes based on orientation of genes
    """
    INTERCHROMOSOMAL = "Interchromosomal"
    OVERLAP_SAME = "Overlapping_Same"
    OVERLAP_CONVERGE = "Overlapping_Converging"
    OVERLAP_DIVERGE = "Overlapping_Diverging"
    OVERLAP_COMPLEX = "Overlapping_Complex"
    READTHROUGH = "Read_Through"
    ADJ_CONVERGE = "Adjacent_Converging"
    ADJ_DIVERGE = "Adjacent_Diverging"
    ADJ_COMPLEX = "Adjacent_Complex"
    INTRACHROMOSOMAL = "Intrachromosomal"
    INTRA_CONVERGE = "Intrachromosomal_Converging"
    INTRA_DIVERGE = "Intrachromsomal_Diverging"
    INTRA_COMPLEX = "Intrachromosomal_Complex"
    UNKNOWN = "Undetermined"

def get_chimera_type(cluster5p, cluster3p, 
                     transcripts5p, transcripts3p, 
                     transcript_dict, genome_tx_trees):
    """
    return tuple containing ChimeraType and distance 
    between 5' and 3' genes 
    """
    chrom5p, start5p, end5p, strand5p = cluster5p.rname, cluster5p.start, cluster5p.end, cluster5p.strand
    chrom3p, start3p, end3p, strand3p = cluster3p.rname, cluster3p.start, cluster3p.end, cluster3p.strand
    # interchromosomal
    if chrom5p != chrom3p:
        return ChimeraTypes.INTERCHROMOSOMAL, -1    
    # strandedness
    same_strand = (strand5p == strand3p)
    # partner orientation
    partners_oriented = ((start5p <= start3p and strand5p == "+") or
                         (start5p > start3p and strand5p == "-")) 
    # genes on same chromosome so check overlap
    is_overlapping = (start5p < end3p) and (start3p < end5p)            
    if is_overlapping:
        if same_strand:
            if partners_oriented:
                return (ChimeraTypes.OVERLAP_SAME, 0)
            else:
                return (ChimeraTypes.OVERLAP_COMPLEX, 0)
        else:
            if partners_oriented:
                return (ChimeraTypes.OVERLAP_CONVERGE, 0)
            else:
                return (ChimeraTypes.OVERLAP_DIVERGE, 0)
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
    my_tx_ids = set([t.tx_id for t in transcripts5p])
    my_tx_ids.update([t.tx_id for t in transcripts3p])
    tx_ids = set()
    for hit in genome_tx_trees[chrom5p].find(between_start, between_end):
        if hit.value in my_tx_ids:
            continue 
        if (hit.start > between_start and
            hit.end < between_end):
            tx_ids.add(hit.value)
    between = []
    between_same_strand = []
    for tx_id in tx_ids:
        t = transcript_dict[tx_id]
        between.append(hit)
        if t.strand == strand5p:
            between_same_strand.append(hit)
    # logic for determining chimera type
    if same_strand:
        if partners_oriented:
            if len(between_same_strand) == 0:
                # genes on same strand, no intervening genes, and 
                # 5' -> 3' partners match gene orientation
                return ChimeraTypes.READTHROUGH, distance
            else:                
                # genes on same strand, partners oriented, but 
                # intervening genes
                return ChimeraTypes.INTRACHROMOSOMAL, distance
        else:
            if len(between) == 0:
                # no intervening genes but partners in opposite orientations
                return (ChimeraTypes.ADJ_COMPLEX, distance)
            else:
                # intervening genes with partners in opposite orientations
                return (ChimeraTypes.INTRA_COMPLEX, distance)
    else:        
        # genes on opposite strands so has to be a complex rearrangement 
        # of some kind
        if len(between) == 0:
            if partners_oriented:
                # 5' -> 3' genomic orientation maintained
                return (ChimeraTypes.ADJ_CONVERGE, distance)
            else:
                return (ChimeraTypes.ADJ_DIVERGE, distance)
        else:
            # intervening genes
            if partners_oriented:
                # 5' -> 3' genomic orientation maintained
                return (ChimeraTypes.INTRA_CONVERGE, distance)
            else:
                return (ChimeraTypes.INTRA_DIVERGE, distance)
