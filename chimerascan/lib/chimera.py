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
    DISCORDANT_SPANNING = 32

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

# chimera types
class ChimeraTypes(object):
    INTERCHROMOSOMAL = "Interchromosomal"
    OVERLAP_CONVERGE = "Overlapping_Converging"
    OVERLAP_DIVERGE = "Overlapping_Diverging"
    OVERLAP_SAME = "Overlapping_Same"
    OVERLAP_COMPLEX = "Overlapping_Complex"
    READTHROUGH = "Read_Through"
    INTRACHROMOSOMAL = "Intrachromosomal"
    ADJ_CONVERGE = "Adjacent_Converging"
    ADJ_DIVERGE = "Adjacent_Diverging"
    ADJ_SAME = "Adjacent_Same"
    ADJ_COMPLEX = "Adjacent_Complex"
    INTRA_COMPLEX = "Intrachromosomal_Complex"
    UNKNOWN = "Undetermined"
    #CHIMERA_INTRA_CONVERGE = "Intrachromosomal_Converging"
    #CHIMERA_INTRA_DIVERGE = "Intrachromsomal_Diverging"
    #CHIMERA_INTRA_SAME = "Intrachromosomal_Same"

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
        return a

    def to_list(self):
        return [self.qname, self.hit_index, self.readnum, self.seq, 
                self.tid, self.pos, self.aend, self.clipstart, 
                self.clipend, int(self.is_reverse), self.numhits, 
                self.mismatches, self.discordant_type, 
                self.orientation]


class ChimeraPartner(object):
    NUM_FIELDS = 12
    
    def __init__(self):
        self.gene_name = None
        self.tx_name = None
        self.start = 0
        self.end = 0
        self.strand = "."
        self.exon_start_num = 0
        self.exon_end_num = 0
        self.inner_dist = 0
        self.mismatches = 0
        self.frac = 0.0
        self.multimap_hist = []
        self.weighted_cov = 0.0
    
    def __repr__(self):
        return ("<%s(gene_name='%s',tx_name='%s')>" % 
                (self.__class__.__name__, self.gene_name, self.tx_name))

    @staticmethod
    def from_list(fields):
        p = ChimeraPartner()
        p.gene_name = fields[0]
        p.tx_name = fields[1]
        p.start = int(fields[2])
        p.end = int(fields[3])
        p.strand = fields[4]
        p.exon_start_num = fields[5]
        p.exon_end_num = fields[6]
        p.inner_dist = fields[7]
        p.mismatches = int(fields[8])
        p.frac = float(fields[9])
        p.multimap_hist = map(int, fields[10].split(","))
        p.weighted_cov = float(fields[11])
        return p
    
    def to_list(self):
        return [self.gene_name, self.tx_name, self.start, self.end, 
                self.strand, self.exon_start_num, self.exon_end_num,
                self.inner_dist, self.mismatches, self.frac, 
                ",".join(map(str,self.multimap_hist)), self.weighted_cov]


class Chimera(object):
    FIELD_DELIM = "|"
    PAIR_DELIM = "||" 
    READ_DELIM = ";"

    def __init__(self):
        self.partner5p = None
        self.partner3p = None
        self.name = None
        self.chimera_type = 0
        self.distance = None
        self.num_encomp_frags = 0
        self.num_spanning_reads = 0
        self.read1_sense_frags = 0
        # junction information
        self.breakpoint_name = None
        self.breakpoint_homology_5p = 0
        self.breakpoint_homology_3p = 0
        # raw read information
        self.encomp_read_pairs = []
        self.spanning_reads = []       

    @staticmethod
    def from_list(fields):
        c = Chimera()
        # chimera partner information
        c.partner5p = ChimeraPartner.from_list(fields[0:ChimeraPartner.NUM_FIELDS])
        c.partner3p = ChimeraPartner.from_list(fields[ChimeraPartner.NUM_FIELDS:2*ChimeraPartner.NUM_FIELDS])
        # chimera information
        FIRSTCOL = 2*ChimeraPartner.NUM_FIELDS
        c.name = fields[FIRSTCOL]
        c.chimera_type = fields[FIRSTCOL+1]
        c.distance = parse_string_none(fields[FIRSTCOL+2])
        if c.distance is not None:
            c.distance = int(c.distance)
        c.num_encomp_frags = int(fields[FIRSTCOL+3])        
        c.num_spanning_reads = int(fields[FIRSTCOL+4])
        c.read1_sense_frags = int(fields[FIRSTCOL+5])
        # breakpoint information
        c.breakpoint_name = parse_string_none(fields[FIRSTCOL+6])
        c.breakpoint_homology_5p = parse_string_none(fields[FIRSTCOL+7])
        if c.breakpoint_homology_5p is not None:
            c.breakpoint_homology_5p = int(c.breakpoint_homology_5p)
        c.breakpoint_homology_3p = parse_string_none(fields[FIRSTCOL+8])
        if c.breakpoint_homology_3p is not None:
            c.breakpoint_homology_3p = int(c.breakpoint_homology_3p)
        # raw encompassing read information
        encomp_reads_field = parse_string_none(fields[FIRSTCOL+9])
        if encomp_reads_field is None:
            c.encomp_read_pairs = []
        for read_pair_fields in encomp_reads_field.split(c.READ_DELIM):
            dreads = []
            for read_fields in read_pair_fields.split(c.PAIR_DELIM):
                dreads.append(DiscordantRead.from_list(read_fields.split(c.FIELD_DELIM)))
            c.encomp_read_pairs.append(dreads)
        # raw spanning read information
        spanning_reads_field = parse_string_none(fields[FIRSTCOL+10])
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
        fields = []
        fields.extend(self.partner5p.to_list())
        fields.extend(self.partner3p.to_list())
        fields.extend([self.name, self.chimera_type, self.distance,
                       self.num_encomp_frags, self.num_spanning_reads,
                       self.read1_sense_frags, 
                       self.breakpoint_name,
                       self.breakpoint_homology_5p, 
                       self.breakpoint_homology_3p])
        encomp_reads = []
        for dreads in self.encomp_read_pairs:
            r5p = self.FIELD_DELIM.join(map(str,dreads[0].to_list()))
            r3p = self.FIELD_DELIM.join(map(str,dreads[1].to_list()))                        
            pair_fields = self.PAIR_DELIM.join([r5p,r3p])
            encomp_reads.append(pair_fields)
        fields.append(self.READ_DELIM.join(encomp_reads))
        spanning_reads = []
        for dread in self.spanning_reads:
            r = self.FIELD_DELIM.join(map(str,dread.to_list()))                        
            spanning_reads.append(r)
        if len(spanning_reads) == 0:
            fields.append("None")
        else:
            fields.append(self.READ_DELIM.join(spanning_reads))
        return fields

    def get_weighted_cov(self):
        """
        weighted coverage is average of weighted coverage of 
        both partners in the chimera
        """
        return 0.5 * (self.partner5p.weighted_cov +
                      self.partner3p.weighted_cov)

