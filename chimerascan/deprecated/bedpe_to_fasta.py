'''
Created on Oct 26, 2010

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
import collections
import os

# local imports
from chimerascan import pysam
from chimerascan.lib import config

def find_homology(seq1, seq2, num_mismatches):
    smallest_len = min(len(seq1), len(seq2))
    mm = 0
    i = 0
    for i in xrange(smallest_len):
        if seq1[i] != seq2[i]:
            mm += 1
            if mm > num_mismatches:
                break
    #print '\t'.join(map(str, [i, seq1, seq2]))
    return i

def bedpe_to_junction_fasta(bed_file, reference_seq_file, read_length,
                            fasta_output_fh, junc_output_fh,
                            num_mismatches=2):
    gene_fasta_prefix = config.GENE_REF_PREFIX
    ref_fa = pysam.Fastafile(reference_seq_file)
    juncs = collections.defaultdict(lambda: [])
    for line in open(bed_file):
        #print line
        fields = line.strip().split('\t')
        ref5p, start5p, end5p = fields[0], int(fields[1]), int(fields[2])
        ref3p, start3p, end3p = fields[3], int(fields[4]), int(fields[5])
        # join end of 5' ref with beginning of 3' ref
        junc_start5p = max(start5p, end5p - read_length + 1)
        junc_end3p = min(end3p, start3p + read_length - 1)
        # fetch sequence
        seq5p = ref_fa.fetch(gene_fasta_prefix + ref5p, junc_start5p, end5p)
        seq3p = ref_fa.fetch(gene_fasta_prefix + ref3p, start3p, junc_end3p)
        seq = seq5p + seq3p
        if len(seq) < (read_length*2) - 2:
            logging.warning("Could not extract sequence of length >%d from BEDPE, only retrieved sequence of (%d,%d) for gene %s" % 
                            ((read_length*2)-2, len(seq5p), len(seq3p), line.strip()))
        # fetch continuation sequence of non-fusion gene
        homolog_end5p = end5p + read_length - 1
        homolog_start3p = max(0, start3p - read_length + 1)
        homolog5p = ref_fa.fetch(gene_fasta_prefix + ref3p, homolog_start3p, start3p)
        homolog3p = ref_fa.fetch(gene_fasta_prefix + ref5p, end5p, homolog_end5p)
        # find homology between 5' gene and 3' gene
        homology_length_5p = find_homology(seq5p, homolog5p, num_mismatches)
        homology_length_3p = find_homology(seq3p, homolog3p, num_mismatches)
        # add sequence to dictionary and group fusion candidates together
        # if they have the same junction sequence
        juncs[seq].append((len(seq5p), homology_length_5p, homology_length_3p, fields))
    # now extract the unique junction sequences
    # and write them to a fasta file
    junc_index = 1    
    for junc_seq,junc_info_list in juncs.iteritems():
        junc_name = "JUNC%07d" % (junc_index) 
        # write to fasta file
        print >>fasta_output_fh, ">%s\n%s" % (junc_name, junc_seq)
        # create entries in junc map file
        for junc_info in junc_info_list:
            left_seq_length, homology_length_5p, homology_length_3p, bedpe_fields = junc_info
            fields = [junc_name, left_seq_length, 
                      homology_length_5p, homology_length_3p]
            fields.extend(bedpe_fields)
            print >>junc_output_fh, '\t'.join(map(str, fields))
        junc_index += 1

def main():
    from optparse import OptionParser
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = OptionParser("usage: %prog [options] <chimeras.bedpe> <out.fasta> <out.juncs>")
    parser.add_option("--rlen", type="int", dest="read_length", default=None)
    parser.add_option("--homology-mismatches", type="int", 
                      dest="num_homology_mismatches", default=2,
                      help="Number of mismatches to tolerate when computing "
                      "homology between gene and its chimeric partner")
    parser.add_option("--index", dest="index_dir") 
    options, args = parser.parse_args()
    bedpe_file = args[0]
    output_fasta_file = args[1]
    output_junc_file = args[2]
    ref_fasta_file = os.path.join(options.index_dir, config.ALIGN_INDEX + ".fa")    
    bedpe_to_junction_fasta(bedpe_file, ref_fasta_file,
                            options.read_length,
                            open(output_fasta_file, "w"),
                            open(output_junc_file, "w"),
                            options.num_homology_mismatches)

if __name__ == '__main__':
    main()

