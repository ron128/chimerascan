'''
Created on Dec 18, 2010

@author: mkiyer
'''
import logging
import itertools

class GeneFeature(object):    
    __slots__ = ('chrom', 'tx_start', 'tx_end', 'tx_name', 'gene_name', 
                 'strand', 'cds_start', 'cds_end', 'exon_count', 'exons') 

    def __str__(self):
        fields = [self.tx_name,
                  self.chrom,
                  self.strand,
                  str(self.tx_start),
                  str(self.tx_end),
                  str(self.cds_start),
                  str(self.cds_end),
                  str(self.exon_count),
                  ','.join(map(str, [e[0] for e in self.exons])) + ',',
                  ','.join(map(str, [e[1] for e in self.exons])) + ',',
                  self.gene_name]
        return '\t'.join(fields)

    @staticmethod
    def from_string(line):        
        if line is None:
            return None
        line = line.strip()
        if line.startswith('#'):
            logging.debug("skipping comment line: %s" % (line))
            return None
        if line.startswith('track'):
            logging.debug("skipping track header line: %s"  % (line))
            return None
        fields = line.split('\t')
        # first six fields are required
        g = GeneFeature()
        g.tx_name = fields[0]
        g.chrom = fields[1]
        g.strand = fields[2]
        g.tx_start = int(fields[3])
        g.tx_end = int(fields[4])
        g.cds_start = int(fields[5])
        g.cds_end = int(fields[6])
        g.exon_count = int(fields[7])
        exon_starts = map(int, fields[8].split(',')[:-1])
        exon_ends = map(int, fields[9].split(',')[:-1])
        g.exons = zip(exon_starts, exon_ends)
        g.gene_name = fields[10]
        return g

    @staticmethod
    def parse(line_iter):
        for line in line_iter:
            if not line:
                continue
            if not line.strip():
                continue            
            if line.startswith("#"):
                continue
            if line.startswith("track"):
                continue
            yield GeneFeature.from_string(line)

class BEDFeature(object):    
    __slots__ = ('chrom', 'tx_start', 'tx_end', 'name', 'score', 'strand',
                 'cds_start', 'cds_end', 'exon_count', 'block_starts', 
                 'block_sizes', 'exons', 'attr_fields')

    def __str__(self):
        fields = [self.chrom,
                  str(self.tx_start),
                  str(self.tx_end),
                  self.name,
                  str(self.score),
                  self.strand,
                  str(self.cds_start),
                  str(self.cds_end),
                  '0',
                  str(self.exon_count),
                  ','.join(map(str, self.block_sizes)) + ',',
                  ','.join(map(str, self.block_starts)) + ',']
        return '\t'.join(fields)

    @staticmethod
    def from_string(line):        
        if line is None:
            return None
        line = line.strip()
        if line.startswith('#'):
            logging.debug("skipping comment line: %s" % (line))
            return None
        if line.startswith('track'):
            logging.debug("skipping track header line: %s"  % (line))
            return None
        fields = line.split('\t')
        # first six fields are required
        g = BEDFeature()
        g.chrom = fields[0]
        g.tx_start = int(fields[1])
        g.tx_end = int(fields[2])
        g.name = fields[3]
        if len(fields) <= 4:
            g.score = 0
            g.strand = '.'
        else:
            g.score = fields[4]
            g.strand = fields[5]        
        if len(fields) <= 6:
            g.cds_start = g.tx_start
            g.cds_end = g.tx_end
            g.exon_count = 1
            g.exons = [(g.tx_start, g.tx_end)]
        else:
            g.cds_start = int(fields[6])
            g.cds_end = int(fields[7])
            g.exon_count = int(fields[9])
            g.block_sizes = map(int, fields[10].split(',')[:-1])
            g.block_starts = map(int, fields[11].split(',')[:-1])            
            g.exons = []
            for start, size in itertools.izip(g.block_starts, g.block_sizes):
                g.exons.append((g.tx_start + start, g.tx_start + start + size))
        if len(fields) <= 12:
            g.attr_fields = []
        else:
            g.attr_fields = fields[12:]
        return g

    @staticmethod
    def parse(line_iter):
        for line in line_iter:
            if not line:
                continue
            if not line.strip():
                continue            
            if line.startswith("#"):
                continue
            if line.startswith("track"):
                continue
            yield BEDFeature.from_string(line)

#class BEDGene():
#    pass
#
#def parse_bed12_line(line):
#    if line is None:
#        return None
#    line = line.strip()
#    if line.startswith('#'):
#        return None
#    if line.startswith('track'):
#        return None
#    fields = line.split('\t')
#    # first six fields are required
#    g = BEDGene()
#    g.chrom = fields[0]
#    g.tx_start = int(fields[1])
#    g.tx_end = int(fields[2])
#    g.name = fields[3]
#    g.score = fields[4]
#    g.strand = fields[5]        
#    g.cds_start = int(fields[6])
#    g.cds_end = int(fields[7])
#    g.exon_count = int(fields[9])
#    block_sizes = map(int, fields[10].split(',')[:-1])
#    block_starts = map(int, fields[11].split(',')[:-1])        
#    g.exon_starts = [(g.tx_start + start) for start in block_starts]        
#    g.exon_ends = [(start + size) for start, size in zip(g.exon_starts, block_sizes)]
#    g.exons = zip(g.exon_starts, g.exon_ends)
#    g.introns = zip(g.exon_ends, g.exon_starts[1:])        
#    return g
#
#def parse_bed12_file(line_iter):
#    '''
#    parse a gene bed file
#    '''
#    for line in line_iter:
#        g = parse_bed12_line(line)
#        if g is None:
#            continue
#        yield g