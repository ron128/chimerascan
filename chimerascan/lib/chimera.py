'''
Created on Jun 3, 2011

@author: mkiyer
'''
from base import parse_string_none
from sam import get_clipped_interval
from stats import hist, scoreatpercentile
import collections

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


class ChimeraPartner(object):
    NUM_FIELDS = 11

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
        p.inner_dist = int(fields[7])
        p.mismatches = int(fields[8])
        p.multimap_hist = map(int, fields[9].split(","))
        p.weighted_cov = float(fields[10])
        return p
    
    def to_list(self):
        return [self.gene_name, self.tx_name, self.start, self.end, 
                self.strand, self.exon_start_num, self.exon_end_num,
                self.inner_dist, self.mismatches, 
                ','.join(map(str, self.multimap_hist)),
                self.weighted_cov]

    @staticmethod
    def from_discordant_reads(dreads, g, 
                              trim_bp=EXON_JUNCTION_TRIM_BP):
        p = ChimeraPartner()
        p.dreads = dreads
        # fix gene names with spaces    
        p.gene_name = '_'.join(g.gene_name.split())
        p.tx_name = g.tx_name
        p.strand = g.strand
        # gather statistics from all the discordant reads that align 
        # to this gene. these statistics will be used later to choose 
        # the most probable chimeras
        starts = []
        ends = []
        multimaps = []
        mismatches = 0
        num5p, num3p = 0,0
        for r in dreads:
            starts.append(r.pos)
            ends.append(r.aend)
            multimaps.append(r.numhits)
            mismatches += r.mismatches
            if r.orientation == OrientationTags.FIVEPRIME:
                num5p += 1
            elif r.orientation == OrientationTags.THREEPRIME:
                num3p += 1
        starts = sorted(starts)
        ends = sorted(ends)
        # get exons corresponding to the start/end, and trim to account for the
        # occasional presence of unmapped "splash" around exon boundaries
        # TODO: the correct way to account for this is to clip reads that have 
        # mismatches to the edge of the nearest exon, accounting for any homology
        # between genomic DNA and the transcript
        firststart, lastend = starts[0], ends[-1]
        trimstart = min(firststart + trim_bp, lastend)
        trimend = max(lastend - trim_bp, trimstart + 1)
        # translate transcript position to exon number
        firstexon_num, firstexon_start, firstexon_end = g.get_exon_interval(trimstart)
        lastexon_num, lastexon_start, lastexon_end = g.get_exon_interval(trimend)
        # set start/end position of chimera partner
        if (num5p > 0) and (num3p > 0):
            raise ValueError("Both 5' and 3' reads found in ChimeraPartner")
        elif num5p > 0:
            p.start = 0
            p.end = lastexon_end
        else:
            tx_length = sum(end - start for start, end in g.exons)
            p.start = firstexon_start
            p.end = tx_length
        p.exon_start_num = firstexon_num
        p.exon_end_num = lastexon_num
        # estimate inner distance between read pairs based on the
        # distribution of reads on this chimera.  inner distance plus
        # twice the segment length should equal the total fragment length
        # TODO: here, we use 95% of start/end positions to account for
        # spuriously large fragments
        p.inner_dist = int(scoreatpercentile(ends, 0.95) - 
                           scoreatpercentile(starts, 0.05))
        # total number of mismatches in all reads
        p.mismatches = mismatches
        # get a histogram of the multimapping profile of this chimera
        p.multimap_hist = hist(multimaps, MULTIMAP_BINS)
        # compute the weighted coverage (allocating equal weight to each
        # mapping of ambiguous reads)
        p.weighted_cov = sum((1.0/x) for x in multimaps) 
        return p


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
        # junction information
        self.breakpoint_name = None
        self.breakpoint_homology_5p = 0
        self.breakpoint_homology_3p = 0
        # raw read information
        self.encomp_read_pairs = []
        self.spanning_reads = []       

    def to_list(self):
        fields = []
        fields.extend(self.partner5p.to_list())
        fields.extend(self.partner3p.to_list())
        fields.extend([self.name, 
                       self.chimera_type, 
                       self.distance,
                       self.breakpoint_name,
                       self.breakpoint_homology_5p, 
                       self.breakpoint_homology_3p])
        # encompassing read pairs
        encomp_reads = []
        for dreads in self.encomp_read_pairs:
            r5p = self.FIELD_DELIM.join(map(str,dreads[0].to_list()))
            r3p = self.FIELD_DELIM.join(map(str,dreads[1].to_list()))                        
            pair_fields = self.PAIR_DELIM.join([r5p,r3p])
            encomp_reads.append(pair_fields)
        fields.append(self.READ_DELIM.join(encomp_reads))
        # spanning single reads
        if len(self.spanning_reads) == 0:
            fields.append("None")
        else:
            spanning_reads = []
            for dread in self.spanning_reads:
                r = self.FIELD_DELIM.join(map(str,dread.to_list()))                        
                spanning_reads.append(r)
            fields.append(self.READ_DELIM.join(spanning_reads))
        return fields

    @staticmethod
    def from_list(fields):
        c = Chimera()
        # 5' and 3' partner information
        c.partner5p = ChimeraPartner.from_list(fields[0:ChimeraPartner.NUM_FIELDS])
        FIRSTCOL = ChimeraPartner.NUM_FIELDS
        c.partner3p = ChimeraPartner.from_list(fields[FIRSTCOL:FIRSTCOL+ChimeraPartner.NUM_FIELDS])
        FIRSTCOL = FIRSTCOL+ChimeraPartner.NUM_FIELDS
        # chimera information
        c.name = fields[FIRSTCOL+0]
        c.chimera_type = fields[FIRSTCOL+1]
        c.distance = parse_string_none(fields[FIRSTCOL+2])
        if c.distance is not None:
            c.distance = int(c.distance)
        # breakpoint information
        c.breakpoint_name = parse_string_none(fields[FIRSTCOL+3])
        c.breakpoint_homology_5p = parse_string_none(fields[FIRSTCOL+4])
        if c.breakpoint_homology_5p is not None:
            c.breakpoint_homology_5p = int(c.breakpoint_homology_5p)
        c.breakpoint_homology_3p = parse_string_none(fields[FIRSTCOL+5])
        if c.breakpoint_homology_3p is not None:
            c.breakpoint_homology_3p = int(c.breakpoint_homology_3p)
        # raw encompassing read information
        encomp_reads_field = parse_string_none(fields[FIRSTCOL+6])
        if encomp_reads_field is None:
            c.encomp_read_pairs = []
        for read_pair_fields in encomp_reads_field.split(c.READ_DELIM):
            dreads = []
            for read_fields in read_pair_fields.split(c.PAIR_DELIM):
                dreads.append(DiscordantRead.from_list(read_fields.split(c.FIELD_DELIM)))
            c.encomp_read_pairs.append(dreads)
        # raw spanning read information
        spanning_reads_field = parse_string_none(fields[FIRSTCOL+7])
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

    def get_weighted_cov(self):
        """
        weighted coverage is the number of reads supporting the
        chimera divided by the number of alignments of the reads
        such that multimapping reads will be assigned a fractional
        weight
        """
        cov = 0.0  
        qnames = set()
        for dpair in self.encomp_read_pairs:
            cov += 2.0 / (dpair[0].numhits + dpair[1].numhits)
            qnames.add(dpair[0].qname)
        for dr in self.spanning_reads:
            if dr.qname not in qnames:
                cov += 1.0 / dr.numhits
        return cov

    def get_weighted_unique_frags(self):
        """
        weighted coverage is the number of unique alignment positions 
        supporting the chimera divided by the number of alignments 
        of the reads such that multimapping reads will be assigned 
        a fractional weight
        """
        encomp_pos_dict = collections.defaultdict(lambda: 1e6)
        qnames = set()
        # add encompassing reads
        for pair in self.encomp_read_pairs:
            key = (pair[0].pos, pair[1].pos)
            weighted_numhits = (pair[0].numhits + pair[1].numhits)/2.0
            encomp_pos_dict[key] = min(encomp_pos_dict[key], weighted_numhits)
            qnames.add(pair[0].qname)
        # add spanning reads
        # TODO: we don't have a way to compute the "true" number of multimapping
        # hits for spanning reads, since breakpoints may be almost the same but
        # have slightly different sequences.  the best thing to do would be to
        # translate spanning reads to their relative position on the transcript
        # and then keep track of unique positions in transcript space
        spanning_pos_set = set()
        for dr in self.spanning_reads:
            if dr.qname not in qnames:
                spanning_pos_set.add(dr.pos)
        # compute weighted unique frags
        wtfrags = sum(1.0/numhits for numhits in encomp_pos_dict.itervalues())
        wtfrags += len(spanning_pos_set)
        return wtfrags

    def get_num_frags(self):
        qnames = set()
        for pair in self.encomp_read_pairs:
            qnames.add(pair[0].qname)
        for dr in self.spanning_reads:
            qnames.add(dr.qname)
        return len(qnames)

    def get_num_spanning_frags(self):
        qnames = set()
        for dr in self.spanning_reads:
            qnames.add(dr.qname)
        return len(qnames)  

    def get_num_unique_positions(self):
        """
        calculates total number of unique read alignment
        positions supporting chimera
        """
        # find all unique alignment positions and read names
        encomp_pos = set()
        qnames = set()
        for pair in self.encomp_read_pairs:
            encomp_pos.add((pair[0].pos, pair[1].pos))
            qnames.add(pair[0].qname)
        # add spanning reads
        spanning_pos = set()
        for dr in self.spanning_reads:
            if dr.qname not in qnames:
                spanning_pos.add(dr.pos)
        return len(encomp_pos) + len(spanning_pos)

    def get_num_unique_spanning_positions(self):
        pos = set()
        for dr in self.spanning_reads:
            pos.add(dr.pos)
        return len(pos)

    def get_read1_sense_frags(self):
        """
        returns the number of fragments where read1 aligns
        in the 'sense' direction relative to the transcript
        """
        return sum(int(r5p.readnum == 0) for r5p,r3p in self.encomp_read_pairs)
