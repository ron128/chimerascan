'''
Created on Jan 31, 2011

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
import unittest

features = ['\t'.join(["uc001aac.3", "chr1","-","14362","29370","14362","14362","11",
                       "14362,14969,15795,16606,16857,17258,17605,17914,18267,24737,29320,",
                       "14829,15038,15947,16765,17055,17368,17742,18061,18369,24891,29370,",
                       "DKFZp434K1323"]),
            '\t'.join(["uc001abw.1", "chr1", "+", "861120", "879961", "861321", "879533", "14",
                       "861120,861301,865534,866418,871151,874419,874654,876523,877515,877789,877938,878632,879077,879287,",
                       "861180,861393,865716,866469,871276,874509,874840,876686,877631,877868,878438,878757,879188,879961,",
                       "SAMD11"])]

from ..gene_to_genome2 import build_gene_to_genome_map, gene_to_genome_pos
from ..feature import GeneFeature

class TestGeneToGenome(unittest.TestCase):

    def test_build(self):
        ggmap = build_gene_to_genome_map(iter(features))        
        for f in features:
            fields = f.strip().split('\t')
            self.assertTrue(fields[0] in ggmap)
            name = fields[0]
            chrom = fields[1]
            # check chrom
            self.assertTrue(ggmap[name][0] == chrom)
            # check strand
            strand = 1 if fields[2] == "-" else 0
            self.assertTrue(ggmap[name][1] == strand)
            # check first/last exon
            start, end = map(int, fields[3:5])
            if strand:
                self.assertTrue(ggmap[name][2][0][1] == end)
                self.assertTrue(ggmap[name][2][-1][0] == start)
            else:
                self.assertTrue(ggmap[name][2][0][0] == start)
                self.assertTrue(ggmap[name][2][-1][1] == end)

    def test_neg_strand(self):
        genes = list(GeneFeature.parse(iter(features)))                
        ggmap = build_gene_to_genome_map(iter(features))        
        chrom,strand,pos = gene_to_genome_pos("uc001aac.3", 0, ggmap)
        self.assertTrue(pos == genes[0].tx_end - 1)
        chrom,strand,pos = gene_to_genome_pos("uc001aac.3", 49, ggmap)
        self.assertTrue(pos == genes[0].tx_end - 50)
        chrom,strand,pos = gene_to_genome_pos("uc001aac.3", 50, ggmap)
        self.assertTrue(pos == genes[0].exons[-2][1] - 1)        
        tx_length = sum(e - s for s,e in genes[0].exons)
        last_pos = tx_length - 1
        chrom,strand,pos = gene_to_genome_pos("uc001aac.3", last_pos, ggmap)
        self.assertTrue(pos == genes[0].tx_start)

    def test_pos_strand(self):
        genes = list(GeneFeature.parse(iter(features)))                
        tx_length = sum(e - s for s,e in genes[1].exons)
        ggmap = build_gene_to_genome_map(iter(features))        
        chrom,strand,pos = gene_to_genome_pos("uc001abw.1", 0, ggmap)
        self.assertTrue(pos == genes[1].tx_start)
        chrom,strand,pos = gene_to_genome_pos("uc001abw.1", 59, ggmap)
        self.assertTrue(pos == genes[1].tx_start + 59)
        chrom,strand,pos = gene_to_genome_pos("uc001abw.1", 60, ggmap)
        self.assertTrue(pos == genes[1].exons[1][0])
        chrom,strand,pos = gene_to_genome_pos("uc001abw.1", 151, ggmap)
        self.assertTrue(pos == genes[1].exons[1][0] + 91)
        chrom,strand,pos = gene_to_genome_pos("uc001abw.1", 152, ggmap)
        self.assertTrue(pos == genes[1].exons[2][0])
        chrom,strand,pos = gene_to_genome_pos("uc001abw.1", tx_length-1, ggmap)
        self.assertTrue(pos == genes[1].tx_end - 1)


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()