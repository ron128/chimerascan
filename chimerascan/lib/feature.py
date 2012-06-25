'''
Created on Dec 18, 2010

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
import itertools
import collections
import operator
import gtf

class TranscriptFeature(object):    
    __slots__ = ('chrom', 'tx_start', 'tx_end', 'strand', 'exon_count',  
                 'exons', 'tx_id', 'cluster_id', 
                 'gene_biotype', 
                 'tx_names',
                 'gene_names', 
                 'annotation_sources')

    def __init__(self):
        self.exons = tuple()
        self.tx_id = -1
        self.cluster_id = -1
        self.gene_biotype = "na"
        self.tx_names = []
        self.gene_names = []
        self.annotation_sources = []
        
    def __str__(self):
        fields = [self.chrom,
                  str(self.tx_start),
                  str(self.tx_end),
                  str(self.tx_id),
                  str(self.cluster_id),
                  self.strand,
                  str(self.exon_count),
                  ','.join(map(str, [e[0] for e in self.exons])) + ',',
                  ','.join(map(str, [e[1] for e in self.exons])) + ',',
                  self.gene_biotype,
                  ','.join(self.tx_names) + ',',
                  ','.join(self.gene_names) + ',',
                  ','.join(self.annotation_sources) + ',']
        return '\t'.join(fields)

    @property
    def introns(self):
        """get tuple of transcript introns"""
        return tuple((self.exons[i-1][1],self.exons[i][0]) 
                     for i in xrange(1,len(self.exons)))
    
    @staticmethod
    def from_string(line):
        if not line:
            return None
        line = line.strip()
        if not line:
            return None
        fields = line.split('\t')
        g = TranscriptFeature()
        g.chrom = fields[0]
        g.tx_start = int(fields[1])
        g.tx_end = int(fields[2])
        g.tx_id = int(fields[3])
        g.cluster_id = int(fields[4])
        g.strand = fields[5]
        g.exon_count = int(fields[6])
        exon_starts = map(int, fields[7].split(',')[:-1])
        exon_ends = map(int, fields[8].split(',')[:-1])
        g.exons = zip(exon_starts, exon_ends)
        g.gene_biotype = fields[9]
        g.tx_names = fields[10].split(',')[:-1]
        g.gene_names = fields[11].split(',')[:-1]
        g.annotation_sources = fields[12].split(',')[:-1]       
        return g
    
    @staticmethod
    def parse(line_iter):
        for line in line_iter:
            g = TranscriptFeature.from_string(line)
            if g is None:
                continue
            yield g

    @staticmethod
    def from_genepred(line_iter):
        cur_transcript_id = 1
        for line in line_iter:
            if not line:
                continue
            if not line.strip():
                continue        
            if not line:
                continue    
            if line.startswith("#"):
                continue
            if line.startswith("track"):
                continue
            fields = line.split('\t')
            # first six fields are required
            g = TranscriptFeature()
            g.tx_names = [fields[0]]
            g.chrom = fields[1]
            g.strand = fields[2]
            g.tx_start = int(fields[3])
            g.tx_end = int(fields[4])
            g.exon_count = int(fields[7])
            exon_starts = map(int, fields[8].split(',')[:-1])
            exon_ends = map(int, fields[9].split(',')[:-1])
            g.exons = zip(exon_starts, exon_ends)
            g.gene_names = [fields[10]]
            g.gene_biotype = "na"
            g.tx_id = cur_transcript_id
            cur_transcript_id += 1
            yield g

    @staticmethod
    def from_bed(line_iter):
        cur_transcript_id = 1
        for line in line_iter:
            if not line:
                continue
            if not line.strip():
                continue        
            if not line:
                continue    
            if line.startswith("#"):
                continue
            if line.startswith("track"):
                continue
            fields = line.split('\t')
            # first six fields are required
            g = TranscriptFeature()
            g.chrom = fields[0]
            g.tx_start = int(fields[1])
            g.tx_end = int(fields[2])
            g.tx_names = [fields[3]]
            g.strand = fields[5]        
            g.exon_count = int(fields[9])
            g.block_sizes = map(int, fields[10].split(',')[:-1])
            g.block_starts = map(int, fields[11].split(',')[:-1])            
            g.exons = []
            for start, size in itertools.izip(g.block_starts, g.block_sizes):
                g.exons.append((g.tx_start + start, g.tx_start + start + size))
            g.gene_biotype = "na"
            g.gene_names = ["na"]
            g.tx_id = cur_transcript_id
            cur_transcript_id += 1
            yield g

    @staticmethod
    def from_gtf(line_iter, source=None):
        chrom_exon_features = collections.defaultdict(lambda: collections.defaultdict(lambda: []))
        for feature in gtf.GTFFeature.parse(line_iter):
            if not feature.feature_type == "exon":
                continue
            if feature.feature_type == "exon":
                transcript_id = feature.attrs["transcript_id"]
                chrom_exon_features[feature.seqid][transcript_id].append(feature)               
        transcripts = []
        cur_transcript_id = 1
        for chrom in sorted(chrom_exon_features):
            exon_features = chrom_exon_features[chrom].values()
            exon_features.sort(key=lambda exon_list: min(x.start for x in exon_list))
            for exons in exon_features:
                # sort exons
                exons.sort(key=operator.attrgetter('start'))
                # build gene feature object
                g = TranscriptFeature()
                g.chrom = exons[0].seqid
                g.tx_start = exons[0].start
                g.tx_end = exons[-1].end
                g.strand = exons[0].strand
                g.exon_count = len(exons)
                g.exons = [(x.start, x.end) for x in exons]
                g.tx_names = [exons[0].attrs['transcript_id']]
                g.gene_names = []
                if 'gene_biotype' in exons[0].attrs:
                    g.gene_biotype = exons[0].attrs['gene_biotype']
                else:
                    g.gene_biotype = 'na'
                if 'gene_name' in exons[0].attrs:
                    g.gene_names.append(exons[0].attrs['gene_name'])
                else:
                    g.gene_names.append('na')
                if source is not None:
                    g.annotation_sources.append(source)
                else:
                    g.annotation_sources.append(exons[0].source)
                g.tx_id = cur_transcript_id 
                cur_transcript_id += 1
                transcripts.append(g)
        return transcripts

