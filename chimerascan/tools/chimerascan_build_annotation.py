#!/usr/bin/env python
'''
Created on May 25, 2012

@author: mkiyer

chimerascan: chimeric transcript discovery using RNA-seq

Copyright (C) 2011-2012 Matthew Iyer

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
import argparse
import sys
import logging
import collections
import operator

from chimerascan.lib import config
from chimerascan.lib.feature import TranscriptFeature
from chimerascan.lib.transcriptome import cluster_transcripts

def build_transcriptome_annotation(gtf_files, output_file):
    # read gtf files and store transcripts
    transcripts = []
    for itm in gtf_files:
        fields = itm.split(",")
        filename,source = fields[0],None
        if len(fields) > 1:
            source = fields[1]
        logging.info("Reading gene features from %s (source=%s)" % (filename, source))
        for t in TranscriptFeature.from_gtf(open(filename), source=source):
            transcripts.append(t)
    logging.debug("\tRead %d annotations from %d files" % (len(transcripts), len(gtf_files)))
    # cluster transcripts by chromosome/strand/position
    logging.info("Determining transcript clusters")
    cur_transcript_id = 1
    cur_cluster_id = 1
    chrom_transcript_clusters = collections.defaultdict(lambda: collections.defaultdict(lambda: []))
    for cluster in cluster_transcripts(transcripts):
        for t in cluster:
            t.tx_id = cur_transcript_id
            t.cluster_id = cur_cluster_id
            chrom_transcript_clusters[t.chrom][(t.introns,t.cluster_id)].append(t)
            cur_transcript_id += 1
        cur_cluster_id += 1
    logging.info("Found %d transcript clusters" % (cur_cluster_id))
    # merge genes in transcript clusters
    logging.info("Merging transcripts")
    outfh = open(output_file, "w")
    cur_transcript_id = 1
    for chrom in sorted(chrom_transcript_clusters):
        transcript_clusters = chrom_transcript_clusters[chrom]
        new_transcripts = []
        for cluster in transcript_clusters.itervalues():
            t = TranscriptFeature()
            t.chrom = chrom
            t.tx_start = min(x.tx_start for x in cluster)
            t.tx_end = max(x.tx_end for x in cluster)
            t.cluster_id = cluster[0].cluster_id 
            t.strand = cluster[0].strand
            t.exon_count = cluster[0].exon_count
            t.exons = list(cluster[0].exons)
            t.exons[0] = (t.tx_start, t.exons[0][1])
            t.exons[-1] = (t.exons[-1][0], t.tx_end)
            t.gene_biotype = "na"
            for x in cluster:
                if x.gene_biotype != "na":
                    t.gene_biotype = x.gene_biotype
                t.tx_names.extend(x.tx_names)
                t.gene_names.extend(x.gene_names)
                t.annotation_sources.extend(x.annotation_sources)
            new_transcripts.append(t)
        new_transcripts.sort(key=operator.attrgetter("tx_start"))
        for t in new_transcripts:
            t.tx_id = cur_transcript_id
            cur_transcript_id += 1
            print >>outfh, str(t)
    outfh.close()
    logging.info("Wrote gene annotation file")

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser(description="Merges multiple gene "
                                     "annotation files into a single "
                                     "transcriptome annotation")
    parser.add_argument("--gtf", dest="gtf_files", action="append")
    parser.add_argument("output_file")
    args = parser.parse_args()
    return build_transcriptome_annotation(args.gtf_files, args.output_file)

if __name__ == '__main__':
    sys.exit(main())
