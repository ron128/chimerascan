'''
Created on Jun 3, 2011

@author: mkiyer
'''
from base import parse_string_none
from sam import get_clipped_interval

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

ORIENTATION_TAG_NAME = "XD"
class OrientationTags(object):
    NONE = 0
    FIVEPRIME = 1
    THREEPRIME = 2

def cmp_orientation(a,b):
    if (a == OrientationTags.NONE) or (b == OrientationTags.NONE):
        return True
    return (a != b)

# constants
MULTIMAP_BINS = (1,2,4,8,16,32,64,128)
CHIMERA_SEP = "|"
# amount of trimming to use to stop reads from overlapping 
# exon boundaries and going into intronic space
EXON_JUNCTION_TRIM_BP = 10

# chimera types
class ChimeraTypes(object):
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

class DiscordantRead(object):
    """
    stores read alignment information needed to nominate 
    chimeric transcripts

    (this is a subset of what is kept in SAM file)
    """
    def __init__(self):
        self.qname = ""
        self.hit_index = -1
        self.readnum = -1
        self.seq = ""
        self.tid = -1
        self.pos = -1
        self.aend = -1
        self.clipstart = -1
        self.clipend = -1
        self.is_reverse = False
        self.numhits = 0
        self.mismatches = 0
        self.discordant_type = 0
        self.orientation = 0
        self.is_spanning = False

    @staticmethod
    def from_read(r):
        a = DiscordantRead()
        a.qname = r.qname
        a.hit_index = r.opt('HI')
        a.readnum = 1 if r.is_read2 else 0
        a.seq = r.seq
        a.tid = r.rname
        a.pos = r.pos
        a.aend = r.aend
        a.clipstart, a.clipend = get_clipped_interval(r)
        a.is_reverse = r.is_reverse
        a.numhits = r.opt('NH')
        a.mismatches = r.opt('NM')
        a.discordant_type = r.opt(DISCORDANT_TAG_NAME)
        a.orientation = r.opt(ORIENTATION_TAG_NAME)
        a.is_spanning = False
        return a

    @staticmethod
    def from_list(fields):
        a = DiscordantRead()
        a.qname = fields[0]
        a.hit_index = int(fields[1])
        a.readnum = int(fields[2])
        a.seq = fields[3]
        a.tid = int(fields[4])
        a.pos = int(fields[5])
        a.aend = int(fields[6])
        a.clipstart = int(fields[7])
        a.clipend = int(fields[8])
        a.is_reverse = True if int(fields[9]) == 1 else False
        a.numhits = int(fields[10])
        a.mismatches = int(fields[11])
        a.discordant_type = int(fields[12])
        a.orientation = int(fields[13])
        a.is_spanning = True if int(fields[14]) == 1 else False
        return a

    def to_list(self):
        return [self.qname, self.hit_index, self.readnum, self.seq, 
                self.tid, self.pos, self.aend, self.clipstart, 
                self.clipend, int(self.is_reverse), self.numhits, 
                self.mismatches, self.discordant_type, 
                self.orientation, int(self.is_spanning)]


def frags_to_encomp_string(frags):
    if len(frags) == 0:
        return "None"
    # encompassing read pairs
    encomp_frags = []
    for frag in frags:
        r5p = Chimera.FIELD_DELIM.join(map(str,frag[0].to_list()))
        r3p = Chimera.FIELD_DELIM.join(map(str,frag[1].to_list()))                        
        pair_fields = Chimera.PAIR_DELIM.join([r5p,r3p])
        encomp_frags.append(pair_fields)
    return Chimera.READ_DELIM.join(encomp_frags)

def get_chimera_type(fiveprime_gene, threeprime_gene, gene_trees):
    """
    return tuple containing ChimeraType and distance 
    between 5' and 3' genes 
    """
    # get gene information
    chrom5p, start5p, end5p, strand5p = fiveprime_gene.chrom, fiveprime_gene.tx_start, fiveprime_gene.tx_end, fiveprime_gene.strand
    chrom3p, start3p, end3p, strand3p = threeprime_gene.chrom, threeprime_gene.tx_start, threeprime_gene.tx_end, threeprime_gene.strand
    # interchromosomal
    if chrom5p != chrom3p:
        return ChimeraTypes.INTERCHROMOSOMAL, None
    # gene orientation
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
    genes_between = []
    genes_between_same_strand = []
    for hit in gene_trees[chrom5p].find(between_start,
                                        between_end):
        if (hit.start > between_start and
            hit.end < between_end):             
            genes_between.append(hit)
            if hit.strand == strand5p:
                genes_between_same_strand.append(hit)
    # logic for determining chimera type
    if same_strand:
        if partners_oriented:
            if len(genes_between_same_strand) == 0:
                # genes on same strand, no intervening genes, and 
                # 5' -> 3' partners match gene orientation
                return ChimeraTypes.READTHROUGH, distance
            else:                
                # genes on same strand, partners oriented, but 
                # intervening genes
                return ChimeraTypes.INTRACHROMOSOMAL, distance
        else:
            if len(genes_between) == 0:
                # no intervening genes but partners in opposite orientations
                return (ChimeraTypes.ADJ_COMPLEX, distance)
            else:
                # intervening genes with partners in opposite orientations
                return (ChimeraTypes.INTRA_COMPLEX, distance)
    else:        
        # genes on opposite strands so has to be a complex rearrangement 
        # of some kind
        if len(genes_between) == 0:
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


class Chimera(object):
    FIELD_DELIM = "|"
    PAIR_DELIM = "||" 
    READ_DELIM = ";"
    TX_NAME_3P_FIELD = 3
    NAME_FIELD = 6
    BREAKPOINT_NAME_FIELD = 14

    def __init__(self):
        self.tx_name_5p = None
        self.tx_start_5p = 0
        self.tx_end_5p = 0
        self.tx_name_3p = None
        self.tx_start_3p = 0
        self.tx_end_3p = 0
        self.name = None
        self.score = 0.0
        self.tx_strand_5p = "."
        self.tx_strand_3p = "."
        self.gene_name_5p = None
        self.gene_name_3p = None
        self.exons_5p = None
        self.exons_3p = None
        self.breakpoint_name = None
        self.breakpoint_seq_5p = None
        self.breakpoint_seq_3p = None
        self.homology_left = None
        self.homology_right = None
        self.encomp_frags = []
        self.spanning_reads = []

    @staticmethod
    def from_list(fields):
        c = Chimera()
        c.tx_name_5p = fields[0]
        c.tx_start_5p = int(fields[1])
        c.tx_end_5p = int(fields[2])
        c.tx_name_3p = fields[3]
        c.tx_start_3p = int(fields[4])
        c.tx_end_3p = int(fields[5])
        c.name = fields[6]
        c.score = fields[7]
        c.tx_strand_5p = fields[8]
        c.tx_strand_3p = fields[9]
        c.gene_name_5p = fields[10]
        c.gene_name_3p = fields[11]
        c.exons_5p = map(int, fields[12].split("-"))
        c.exons_3p = map(int, fields[13].split("-"))
        c.breakpoint_name = fields[14]
        c.breakpoint_seq_5p = fields[15]
        c.breakpoint_seq_3p = fields[16]
        c.homology_left = int(fields[17])
        c.homology_right = int(fields[18])
        c.encomp_frags = []
        c.spanning_reads = []
        # raw encompassing read information
        encomp_reads_field = parse_string_none(fields[19])
        if encomp_reads_field is not None:
            for read_pair_fields in encomp_reads_field.split(c.READ_DELIM):
                dreads = []
                for read_fields in read_pair_fields.split(c.PAIR_DELIM):
                    dreads.append(DiscordantRead.from_list(read_fields.split(c.FIELD_DELIM)))
                c.encomp_frags.append(dreads)
        # raw spanning read information
        spanning_reads_field = parse_string_none(fields[20])
        if spanning_reads_field is not None:
            for read_fields in spanning_reads_field.split(c.READ_DELIM):
                c.spanning_reads.append(DiscordantRead.from_list(read_fields.split(c.FIELD_DELIM)))        
        return c

    @staticmethod
    def parse(line_iter):
        for line in line_iter:
            if line.startswith("#"):
                continue            
            fields = line.strip().split('\t')
            yield Chimera.from_list(fields)

    def to_list(self):
        # reads
        if len(self.spanning_reads) == 0:
            span_string = None
        else:
            span_string = Chimera.READ_DELIM.join(Chimera.FIELD_DELIM.join(map(str,r.to_list())) 
                                                  for r in self.spanning_reads)
        return [self.tx_name_5p, self.tx_start_5p, self.tx_end_5p,
                self.tx_name_3p, self.tx_start_3p, self.tx_end_3p,
                self.name, self.score, 
                self.tx_strand_5p, self.tx_strand_3p,
                self.gene_name_5p, self.gene_name_3p,
                "%d-%d" % (self.exons_5p[0], self.exons_5p[1]),
                "%d-%d" % (self.exons_3p[0], self.exons_3p[1]),
                self.breakpoint_name,
                self.breakpoint_seq_5p,
                self.breakpoint_seq_3p,
                self.homology_left,
                self.homology_right,
                frags_to_encomp_string(self.encomp_frags),
                span_string]

    def get_num_unique_positions(self):
        """
        calculates total number of unique read alignment
        positions supporting chimera
        """
        # find all unique alignment positions and read names
        encomp_pos = set()
        qnames = set()
        for pair in self.encomp_frags:
            if pair[0].qname not in qnames:
                qnames.add(pair[0].qname)
                encomp_pos.add((pair[0].pos, pair[1].pos))
        # add spanning reads
        spanning_pos = set()
        for dr in self.spanning_reads:
            if dr.qname not in qnames:
                qnames.add(dr.qname)
                spanning_pos.add(dr.pos)
        return len(encomp_pos) + len(spanning_pos)

    def get_num_frags(self, maxnumhits=0):
        """
        number of unique fragments supporting the 
        chimera (by read name)
        """
        qnames = set()
        for pair in self.encomp_frags:
            if (maxnumhits > 0) and (min(pair[0].numhits, pair[1].numhits) > maxnumhits):
                continue
            qnames.add(pair[0].qname)
        for dr in self.spanning_reads:
            if (maxnumhits > 0) and (dr.numhits > maxnumhits):
                continue
            qnames.add(dr.qname)
        return len(qnames)

    def get_num_spanning_frags(self, maxnumhits=0):
        """
        number of unique spanning fragments supporting the 
        chimera (by read name)
        """
        qnames = set()
        for dpair in self.encomp_frags:
            if (maxnumhits > 0) and (min(dpair[0].numhits, dpair[1].numhits) > maxnumhits):
                continue
            if any(dr.is_spanning for dr in dpair):
                qnames.add(dpair[0].qname)            
        for dr in self.spanning_reads:
            if (maxnumhits > 0) and (dr.numhits > maxnumhits):
                continue
            qnames.add(dr.qname)
        return len(qnames)  

    def get_spanning_reads(self):
        for dpair in self.encomp_frags:
            if dpair[0].is_spanning:
                yield dpair[0]
            if dpair[1].is_spanning:
                yield dpair[1]
        for dr in self.spanning_reads:
            yield dr
