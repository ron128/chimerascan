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
import logging

# local imports
from feature import GeneFeature

def build_gene_to_genome_map(line_iter, rname_prefix=None):
    # create arrays to map genes in bed file to genome 
    rname_prefix = '' if rname_prefix is None else rname_prefix
    gene_genome_map = {}    
    for g in GeneFeature.parse(line_iter):
        rname = rname_prefix + g.tx_name
        strand = 1 if g.strand == '-' else 0 
        exon_vectors = [(start, end) for start, end in g.exons]
        if strand:
            exon_vectors.reverse()
        if rname in gene_genome_map:
            logging.error("Duplicate references %s found in bed file" % (rname))
        gene_genome_map[rname] = (g.chrom, strand, exon_vectors)
    return gene_genome_map

def gene_to_genome_pos(rname, pos, gene_genome_map):    
    '''
    translate gene 'rname' position 'gene_pos' to genomic
    coordinates.  returns a 3-tuple with (chrom, strand, pos)
    '''
    chrom, strand, intervals = gene_genome_map[rname]
    offset = 0
    for start, end, in intervals:
        exon_size = end - start
        if pos < offset + exon_size:            
            if strand:
                return chrom, strand, start + exon_size - (pos - offset) - 1
            else:
                return chrom, strand, start + (pos - offset)
        #print start, end, offset, pos
        offset += exon_size
    return None
