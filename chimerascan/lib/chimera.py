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
        #print [Chimera.FIELD_DELIM.join(map(str,r.to_list())) for r in self.spanning_reads]
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
            encomp_pos.add((pair[0].pos, pair[1].pos))
            qnames.add(pair[0].qname)
        # add spanning reads
        spanning_pos = set()
        for dr in self.spanning_reads:
            if dr.qname not in qnames:
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
        for dr in self.spanning_reads:
            if (maxnumhits > 0) and (dr.numhits > maxnumhits):
                continue
            qnames.add(dr.qname)
        return len(qnames)  


#
#class ChimeraPartner(object):
#    NUM_FIELDS = 11
#
#    def __init__(self):
#        self.gene_name = None
#        self.tx_name = None
#        self.start = 0
#        self.end = 0
#        self.strand = "."
#        self.exon_start_num = 0
#        self.exon_end_num = 0
#        self.inner_dist = 0
#        self.mismatches = 0
#        self.multimap_hist = []
#        self.weighted_cov = 0.0
#    
#    def __repr__(self):
#        return ("<%s(gene_name='%s',tx_name='%s')>" % 
#                (self.__class__.__name__, self.gene_name, self.tx_name))
#
#    @staticmethod
#    def from_list(fields):
#        p = ChimeraPartner()
#        p.gene_name = fields[0]
#        p.tx_name = fields[1]
#        p.start = int(fields[2])
#        p.end = int(fields[3])
#        p.strand = fields[4]
#        p.exon_start_num = fields[5]
#        p.exon_end_num = fields[6]
#        p.inner_dist = int(fields[7])
#        p.mismatches = int(fields[8])
#        p.multimap_hist = map(int, fields[9].split(","))
#        p.weighted_cov = float(fields[10])
#        return p
#    
#    def to_list(self):
#        return [self.gene_name, self.tx_name, self.start, self.end, 
#                self.strand, self.exon_start_num, self.exon_end_num,
#                self.inner_dist, self.mismatches, 
#                ','.join(map(str, self.multimap_hist)),
#                self.weighted_cov]
#
#    @staticmethod
#    def from_discordant_reads(dreads, g, 
#                              trim_bp=EXON_JUNCTION_TRIM_BP):
#        p = ChimeraPartner()
#        p.dreads = dreads
#        # fix gene names with spaces    
#        p.gene_name = '_'.join(g.gene_name.split())
#        p.tx_name = g.tx_name
#        p.strand = g.strand
#        # gather statistics from all the discordant reads that align 
#        # to this gene. these statistics will be used later to choose 
#        # the most probable chimeras
#        starts = []
#        ends = []
#        multimaps = []
#        mismatches = 0
#        num5p, num3p = 0,0
#        for r in dreads:
#            starts.append(r.pos)
#            ends.append(r.aend)
#            multimaps.append(r.numhits)
#            mismatches += r.mismatches
#            if r.orientation == OrientationTags.FIVEPRIME:
#                num5p += 1
#            elif r.orientation == OrientationTags.THREEPRIME:
#                num3p += 1
#        starts = sorted(starts)
#        ends = sorted(ends)
#        # get exons corresponding to the start/end, and trim to account for the
#        # occasional presence of unmapped "splash" around exon boundaries
#        # TODO: the correct way to account for this is to clip reads that have 
#        # mismatches to the edge of the nearest exon, accounting for any homology
#        # between genomic DNA and the transcript
#        firststart, lastend = starts[0], ends[-1]
#        trimstart = min(firststart + trim_bp, lastend)
#        trimend = max(lastend - trim_bp, trimstart + 1)
#        # translate transcript position to exon number
#        firstexon_num, firstexon_start, firstexon_end = g.get_exon_interval(trimstart)
#        lastexon_num, lastexon_start, lastexon_end = g.get_exon_interval(trimend)
#        # set start/end position of chimera partner
#        if (num5p > 0) and (num3p > 0):
#            raise ValueError("Both 5' and 3' reads found in ChimeraPartner")
#        elif num5p > 0:
#            p.start = 0
#            p.end = lastexon_end
#        else:
#            tx_length = sum(end - start for start, end in g.exons)
#            p.start = firstexon_start
#            p.end = tx_length
#        p.exon_start_num = firstexon_num
#        p.exon_end_num = lastexon_num
#        # estimate inner distance between read pairs based on the
#        # distribution of reads on this chimera.  inner distance plus
#        # twice the segment length should equal the total fragment length
#        # TODO: here, we use 95% of start/end positions to account for
#        # spuriously large fragments
#        p.inner_dist = int(scoreatpercentile(ends, 0.95) - 
#                           scoreatpercentile(starts, 0.05))
#        # total number of mismatches in all reads
#        p.mismatches = mismatches
#        # get a histogram of the multimapping profile of this chimera
#        p.multimap_hist = hist(multimaps, MULTIMAP_BINS)
#        # compute the weighted coverage (allocating equal weight to each
#        # mapping of ambiguous reads)
#        p.weighted_cov = sum((1.0/x) for x in multimaps) 
#        return p
#
#
#
#
#class Chimera(object):
#
#    def __init__(self):
#        self.partner5p = None
#        self.partner3p = None
#        self.name = None
#        self.chimera_type = 0
#        self.distance = None
#        # junction information
#        self.breakpoint_name = None
#        self.breakpoint_homology_5p = 0
#        self.breakpoint_homology_3p = 0
#        # raw read information
#        self.encomp_read_pairs = []
#        self.spanning_reads = []       
#
#    def to_list(self):
#        fields = []
#        fields.extend(self.partner5p.to_list())
#        fields.extend(self.partner3p.to_list())
#        fields.extend([self.name, 
#                       self.chimera_type, 
#                       self.distance,
#                       self.breakpoint_name,
#                       self.breakpoint_homology_5p, 
#                       self.breakpoint_homology_3p])
#        # encompassing read pairs
#        encomp_reads = []
#        for dreads in self.encomp_read_pairs:
#            r5p = self.FIELD_DELIM.join(map(str,dreads[0].to_list()))
#            r3p = self.FIELD_DELIM.join(map(str,dreads[1].to_list()))                        
#            pair_fields = self.PAIR_DELIM.join([r5p,r3p])
#            encomp_reads.append(pair_fields)
#        fields.append(self.READ_DELIM.join(encomp_reads))
#        # spanning single reads
#        if len(self.spanning_reads) == 0:
#            fields.append("None")
#        else:
#            spanning_reads = []
#            for dread in self.spanning_reads:
#                r = self.FIELD_DELIM.join(map(str,dread.to_list()))                        
#                spanning_reads.append(r)
#            fields.append(self.READ_DELIM.join(spanning_reads))
#        return fields
#
#    @staticmethod
#    def from_list(fields):
#        c = Chimera()
#        # 5' and 3' partner information
#        c.partner5p = ChimeraPartner.from_list(fields[0:ChimeraPartner.NUM_FIELDS])
#        FIRSTCOL = ChimeraPartner.NUM_FIELDS
#        c.partner3p = ChimeraPartner.from_list(fields[FIRSTCOL:FIRSTCOL+ChimeraPartner.NUM_FIELDS])
#        FIRSTCOL = FIRSTCOL+ChimeraPartner.NUM_FIELDS
#        # chimera information
#        c.name = fields[FIRSTCOL+0]
#        c.chimera_type = fields[FIRSTCOL+1]
#        c.distance = parse_string_none(fields[FIRSTCOL+2])
#        if c.distance is not None:
#            c.distance = int(c.distance)
#        # breakpoint information
#        c.breakpoint_name = parse_string_none(fields[FIRSTCOL+3])
#        c.breakpoint_homology_5p = parse_string_none(fields[FIRSTCOL+4])
#        if c.breakpoint_homology_5p is not None:
#            c.breakpoint_homology_5p = int(c.breakpoint_homology_5p)
#        c.breakpoint_homology_3p = parse_string_none(fields[FIRSTCOL+5])
#        if c.breakpoint_homology_3p is not None:
#            c.breakpoint_homology_3p = int(c.breakpoint_homology_3p)
#        # raw encompassing read information
#        encomp_reads_field = parse_string_none(fields[FIRSTCOL+6])
#        if encomp_reads_field is None:
#            c.encomp_read_pairs = []
#        for read_pair_fields in encomp_reads_field.split(c.READ_DELIM):
#            dreads = []
#            for read_fields in read_pair_fields.split(c.PAIR_DELIM):
#                dreads.append(DiscordantRead.from_list(read_fields.split(c.FIELD_DELIM)))
#            c.encomp_read_pairs.append(dreads)
#        # raw spanning read information
#        spanning_reads_field = parse_string_none(fields[FIRSTCOL+7])
#        if spanning_reads_field is not None:
#            for read_fields in spanning_reads_field.split(c.READ_DELIM):
#                c.spanning_reads.append(DiscordantRead.from_list(read_fields.split(c.FIELD_DELIM)))
#        return c
#
#    @staticmethod
#    def parse(line_iter):
#        for line in line_iter:
#            if line.startswith("#"):
#                continue            
#            fields = line.strip().split('\t')
#            yield Chimera.from_list(fields)
#
#    def get_weighted_cov(self):
#        """
#        weighted coverage is the number of reads supporting the
#        chimera divided by the number of alignments of the reads
#        such that multimapping reads will be assigned a fractional
#        weight
#        """
#        cov = 0.0  
#        qnames = set()
#        for dpair in self.encomp_read_pairs:
#            cov += 2.0 / (dpair[0].numhits + dpair[1].numhits)
#            qnames.add(dpair[0].qname)
#        for dr in self.spanning_reads:
#            if dr.qname not in qnames:
#                cov += 1.0 / dr.numhits
#        return cov
#
#    def get_weighted_unique_frags(self):
#        """
#        weighted coverage is the number of unique alignment positions 
#        supporting the chimera divided by the number of alignments 
#        of the reads such that multimapping reads will be assigned 
#        a fractional weight
#        """
#        encomp_pos_dict = collections.defaultdict(lambda: 1e6)
#        qnames = set()
#        # add encompassing reads
#        for pair in self.encomp_read_pairs:
#            key = (pair[0].pos, pair[1].pos)
#            weighted_numhits = (pair[0].numhits + pair[1].numhits)/2.0
#            encomp_pos_dict[key] = min(encomp_pos_dict[key], weighted_numhits)
#            qnames.add(pair[0].qname)
#        # add spanning reads
#        # TODO: we don't have a way to compute the "true" number of multimapping
#        # hits for spanning reads, since breakpoints may be almost the same but
#        # have slightly different sequences.  the best thing to do would be to
#        # translate spanning reads to their relative position on the transcript
#        # and then keep track of unique positions in transcript space
#        spanning_pos_set = set()
#        for dr in self.spanning_reads:
#            if dr.qname not in qnames:
#                spanning_pos_set.add(dr.pos)
#        # compute weighted unique frags
#        wtfrags = sum(1.0/numhits for numhits in encomp_pos_dict.itervalues())
#        wtfrags += len(spanning_pos_set)
#        return wtfrags
#
#    def get_num_frags(self):
#        qnames = set()
#        for pair in self.encomp_read_pairs:
#            qnames.add(pair[0].qname)
#        for dr in self.spanning_reads:
#            qnames.add(dr.qname)
#        return len(qnames)
#
#    def get_num_spanning_frags(self):
#        qnames = set()
#        for dr in self.spanning_reads:
#            qnames.add(dr.qname)
#        return len(qnames)  
#
#    def get_num_unique_positions(self):
#        """
#        calculates total number of unique read alignment
#        positions supporting chimera
#        """
#        # find all unique alignment positions and read names
#        encomp_pos = set()
#        qnames = set()
#        for pair in self.encomp_read_pairs:
#            encomp_pos.add((pair[0].pos, pair[1].pos))
#            qnames.add(pair[0].qname)
#        # add spanning reads
#        spanning_pos = set()
#        for dr in self.spanning_reads:
#            if dr.qname not in qnames:
#                spanning_pos.add(dr.pos)
#        return len(encomp_pos) + len(spanning_pos)
#
#    def get_num_unique_spanning_positions(self):
#        pos = set()
#        for dr in self.spanning_reads:
#            pos.add(dr.pos)
#        return len(pos)
#
#    def get_read1_sense_frags(self):
#        """
#        returns the number of fragments where read1 aligns
#        in the 'sense' direction relative to the transcript
#        """
#        return sum(int(r5p.readnum == 0) for r5p,r3p in self.encomp_read_pairs)
