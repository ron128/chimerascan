'''
Created on May 18, 2011

@author: mkiyer
'''

class DiscordantFragment(object):    
    def __init__(self):
        self.qname = None
        self.r1_hit_index = None
        self.r2_hit_index = None
        

class Chimera(object):
    def __init__(self):
        self.id = None
        self.chimera_type = 0
        self.distance = None 
        self.encompassing_reads = 0
        self.weighted_cov = 0.0
        self.strand_reads = [0, 0]
        self.multimap_cov_hist = None
        self.qnames = []
        self.seqs = []       
        self.mate5p = ChimeraMate()
        self.mate3p = ChimeraMate()
        self.extra_fields = []


class ChimeraMate(object):
    def __init__(self):
        self.start = 0
        self.end = 0
        self.tx_name = None
        self.gene_name = None
        self.strand = "."
        self.exon_start_num = 0
        self.exon_end_num = 0
        self.isize = 0
        self.frac = 0.0
                
class Chimera(object):
    QNAME_COL = 22
    LAST_COL = 24
    SEQ_FIELD_DELIM = ';'

    def __init__(self):
        self.name = None
        self.chimera_type = 0
        self.distance = None 
        self.encompassing_reads = 0
        self.weighted_cov = 0.0
        self.strand_reads = [0, 0]
        self.multimap_cov_hist = None
        self.qnames = []
        self.seqs = []       
        self.mate5p = ChimeraMate()
        self.mate3p = ChimeraMate()
        self.extra_fields = []
        




class TranscriptAlignment(object):
    def __init__(self, g, exon_num, chrom, start, end, e_start, e_end, 
                 e_start_overhang, e_end_overhang, 
                 tx_start, tx_end):
        self.g = g
        self.exon_num = exon_num
        self.chrom = chrom
        self.start = start
        self.end = end
        self.e_start = e_start
        self.e_end = e_end
        self.e_start_overhang = e_start_overhang
        self.e_end_overhang = e_end_overhang
        self.tx_start = tx_start
        self.tx_end = tx_end

    def __repr__(self):
        return ("<%s(g='%s', exon_num='%d', chrom='%s', start='%d', "
                "end='%d', e_start='%d', e_end='%d', "
                "e_start_overhang='%d', e_end_overhang='%d', "
                "tx_start='%d', tx_end='%d')>" %
                (self.__class__.__name__, self.g, self.exon_num, 
                 self.chrom, self.start, self.end, self.e_start, self.e_end, 
                 self.e_start_overhang, self.e_end_overhang,
                 self.tx_start, self.tx_end))





class DiscordantCluster(object):
    __slots__ = ('rname', 'start', 'end', 'strand', 'pad_start', 'pad_end', 
                 'multimaps', 'seq', 'qual')
    def __init__(self, rname, start, end, strand, pad_start, pad_end, multimaps,
                 seq=None, 
                 qual=None):
        self.rname = rname
        self.start = start
        self.end = end
        self.strand = strand
        self.pad_start = pad_start
        self.pad_end = pad_end
        self.multimaps = multimaps
        self.seq = seq
        self.qual = qual
    def __repr__(self):
        return ("<%s(rname=%s, strand=%s, start=%d, end=%d, pad_start=%d, "
                "pad_end=%d, multimaps=%d, seq=%s, qual=%s)>" %
                (self.__class__.__name__, self.rname, self.strand, self.start, self.end,
                 self.pad_start, self.pad_end, self.multimaps, self.seq, self.qual))
    def to_list(self):
        return [self.rname, self.start, self.end, self.strand, 
                self.pad_start, self.pad_end, self.multimaps,
                self.seq, self.qual]
    @staticmethod
    def from_list(fields):        
        seq = parse_string_none(fields[7])
        qual = parse_string_none(fields[8])
        return DiscordantCluster(fields[0], int(fields[1]), int(fields[2]), 
                                 fields[3], int(fields[4]), int(fields[5]),
                                 int(fields[6]), seq, qual)


class DiscordantFragment(object):
    __slots__ = ('qname', 'is_genome', 'discordant5p', 'discordant3p',
                 'code', 'read1_is_sense', 'clust5p', 'clust3p')
    # columns in the tabular output that contain the references
    # (useful for sorting)
    REF1_COL = 6
    REF2_COL = 15
    
    # fragment type codes
    NA = 0
    NONMAPPING = 1
    CONCORDANT_SINGLE = 2
    DISCORDANT_SINGLE = 3
    DISCORDANT_SINGLE_COMPLEX = 4
    CONCORDANT_PAIRED = 5
    DISCORDANT_PAIRED = 6
    DISCORDANT_PAIRED_COMPLEX = 7
    _discordant_codes = ["NA", 
                         "NONMAPPING",
                         "CONCORDANT_SINGLE",
                         "DISCORDANT_SINGLE",
                         "DISCORDANT_SINGLE_COMPLEX",
                         "CONCORDANT_PAIRED",
                         "DISCORDANT_PAIRED",
                         "DISCORDANT_PAIRED_COMPLEX"]
    # flag bits
    FLAG_DISCORDANT_5P = 0b10
    FLAG_DISCORDANT_3P = 0b1

    def __init__(self, qname, read1_is_sense, clust5p, clust3p,
                 code=0, is_genome=False, discordant5p=False, 
                 discordant3p=False): 
        self.qname = qname
        self.read1_is_sense = read1_is_sense
        self.clust5p = clust5p
        self.clust3p = clust3p
        self.code = code
        self.is_genome = is_genome
        self.discordant5p = discordant5p
        self.discordant3p = discordant3p

    def __repr__(self):
        return ("<%s(qname=%s, read1_is_sense=%s, clust5p=%s, clust3p=%s "
                "is_genome=%s, discordant5p=%s, discordant3p=%s, "
                "code=%d, string_code=%s)>" %
                (self.__class__.__name__, self.qname, self.read1_is_sense,
                 self.clust5p, self.clust3p, 
                 self.is_genome, self.discordant5p, self.discordant3p, 
                 self.code, self._discordant_codes[self.code]))

    @property
    def clust1(self):
        return self.clust5p if self.read1_is_sense else self.clust3p
    @property
    def clust2(self):
        return self.clust3p if self.read1_is_sense else self.clust5p

    def set_flags(self, nclusts1, nclusts2, nclustspe):
        # set 5'/3' discordant flags according to sense/antisense
        # orientation and cluster mappings
        if self.read1_is_sense:
            self.discordant5p = (nclusts1 > 1)
            self.discordant3p = (nclusts2 > 1)
        else:
            self.discordant5p = (nclusts2 > 1)
            self.discordant3p = (nclusts1 > 1)                    
        if nclustspe == 0:
            self.code = self.NONMAPPING
        elif nclustspe == 1:
            if (nclusts1 == 0) or (nclusts2 == 0):
                # one of the reads is concordant and the other is unmapped
                self.code = self.CONCORDANT_SINGLE
            else:
                # both reads are mapped and concordant
                self.code = self.CONCORDANT_PAIRED
        elif nclustspe == 2:
            if (nclusts1 > 1) or (nclusts2 > 1):
                # one of the reads is discordant and the other is
                # unmapped
                if (nclusts1 > 2) or (nclusts2 > 2):        
                    self.code = self.DISCORDANT_SINGLE_COMPLEX
                else:
                    self.code = self.DISCORDANT_SINGLE
            else:
                self.code = self.DISCORDANT_PAIRED
        else:
            self.code = self.DISCORDANT_PAIRED_COMPLEX
        return self

    def _pack_flags(self):
        '''pack flags into a single integer value for reading/writing to file'''
        flags = ((self.discordant5p << self.FLAG_DISCORDANT_5P) |
                 (self.discordant3p << self.FLAG_DISCORDANT_3P))
        return flags
    
    def _unpack_flags(self, val):
        '''unpack flags and set attributes'''
        self.discordant5p = bool(val & self.FLAG_DISCORDANT_5P)
        self.discordant3p = bool(val & self.FLAG_DISCORDANT_3P)

    def to_list(self):
        genome = "GENOME" if self.is_genome else "GENE"
        return ([self.qname, int(self.read1_is_sense), genome,
                 self._discordant_codes[self.code], self._pack_flags()] +
                self.clust5p.to_list() + self.clust3p.to_list())

    @staticmethod
    def from_list(fields):
        qname = fields[0]        
        read1_is_sense = bool(int(fields[1]))
        is_genome = True if fields[2] == "GENOME" else False
        code = DiscordantFragment._discordant_codes.index(fields[3])
        flags = int(fields[4])
        clust5p = DiscordantCluster.from_list(fields[5:14])
        clust3p = DiscordantCluster.from_list(fields[14:23])
        f = DiscordantFragment(qname, read1_is_sense, clust5p, clust3p,
                               code=code, is_genome=is_genome)
        f._unpack_flags(flags)
        return f
    