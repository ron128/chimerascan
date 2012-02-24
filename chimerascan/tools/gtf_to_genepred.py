'''
Created on Feb 6, 2012

@author: mkiyer

chimerascan: chimeric transcript discovery using RNA-seq

Copyright (C) 2012 Matthew Iyer

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
import operator
import os
import sys
from optparse import OptionParser

from chimerascan.lib import gtf

def gtf_to_genepred(gtf_file, genepred_file):
    # group by transcript id
    logging.info("Reading GTF file")
    chrom_exon_features = collections.defaultdict(lambda: collections.defaultdict(lambda: []))
    for feature in gtf.GTFFeature.parse(open(gtf_file)):
        if feature.feature_type == "exon":
            transcript_id = feature.attrs["transcript_id"]
            chrom_exon_features[feature.seqid][transcript_id].append(feature)
    # convert to genepred
    logging.info("Writing GenePred file")
    outfh = open(genepred_file, "w")
    for chrom in sorted(chrom_exon_features):
        logging.debug("Chromosome %s" % (chrom))
        exon_features = chrom_exon_features[chrom].values()
        exon_features.sort(key=lambda exon_list: min(x.start for x in exon_list))
        for exons in exon_features:
            # sort exons
            exons.sort(key=operator.attrgetter('start'))
            chrom = exons[0].seqid
            tx_start = exons[0].start
            tx_end = exons[-1].end
            strand = exons[0].strand
            transcript_id = exons[0].attrs['transcript_id']
            gene_name = exons[0].attrs['gene_name']
            # write genepred fields
            fields = [transcript_id, chrom, strand, str(tx_start), 
                      str(tx_end), str(tx_start), str(tx_start),
                      str(len(exons)),
                      ",".join(map(str,[x.start for x in exons])) + ",",
                      ",".join(map(str,[x.end for x in exons])) + ",",
                      gene_name]
            print >>outfh, "\t".join(fields)
    outfh.close()

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = OptionParser("usage: %prog <input.gtf> <genepred_output.txt>")
    options, args = parser.parse_args()
    # check command line arguments
    if len(args) < 2:
        parser.error("Incorrect number of command line arguments")
    gtf_file = args[0]
    genepred_file = args[1]
    # check that input files exist
    if not os.path.isfile(gtf_file):
        parser.error("GTF file '%s' not found" % (gtf_file))
    gtf_to_genepred(gtf_file, genepred_file)
    return 0

if __name__ == '__main__':
    sys.exit(main())
