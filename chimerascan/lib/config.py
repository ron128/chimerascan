'''
Created on Jan 5, 2011

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
JOB_SUCCESS = 0
JOB_ERROR = 1

LOG_DIR = "log"
TMP_DIR = "tmp"

# chimerascan index definitions
ALIGN_INDEX = 'align_index'
BOWTIE_INDEX_FILE = 'align_index.1.ebwt'
FRAG_SIZE_INDEX = 'frag_size_index'
FRAG_SIZE_INDEX_FILE = 'frag_size_index.1.ebwt'
GENE_REF_PREFIX = 'gene_'
GENE_FEATURE_FILE = "gene_features.txt"
TOPHAT_JUNCS_FILE = "known_juncs.txt"

# chimerascan run output file definitions
RUNCONFIG_XML_FILE = "runconfig.xml"
TOPHAT_DIR = "tophat"
SORTED_FASTQ_FILES = ("rname_sorted_1.fq", "rname_sorted_2.fq")

FRAG_SIZE_BAM_FILE = "frag_size_reads.bam"
FRAG_SIZE_DIST_FILE = "frag_size_dist.txt"
FRAG_SIZE_MAX_SAMPLES = 1e6
FRAG_SIZE_NUM_STDEVS = 3

ALIGNED_READS_BAM_FILE = "aligned_reads.bam"

UNALIGNED_FASTQ_PARAM = "unaligned.fq"
UNALIGNED_FASTQ_FILES = ("unaligned_1.fq", "unaligned_2.fq")
MAXMULTIMAP_FASTQ_PARAM = "maxmulti.fq"
MAXMULTIMAP_FASTQ_FILES = ("maxmulti_1.fq", "maxmulti_2.fq")
BASE_PROCESSORS = 2
MIN_SEGMENT_LENGTH = 20
DISCORDANT_BAM_FILE = "discordant_reads.bam"
DISCORDANT_PAIRED_BAM_FILE = "discordant_reads_paired.bam"
DISCORDANT_GENE_BEDPE_FILE = "discordant_gene_reads.bedpe"
DISCORDANT_GENOME_BEDPE_FILE = "discordant_genome_reads.bedpe"
EXTENDED_DISCORDANT_GENE_BEDPE_FILE = "discordant_gene_reads.extended.bedpe"
SORTED_DISCORDANT_GENE_BEDPE_FILE = "discordant_gene_reads.srt.bedpe"
SPANNING_FASTQ_FILE = "putative_spanning_reads.fq"
ENCOMPASSING_CHIMERA_BEDPE_FILE = "encompassing_chimeras.bedpe"
FILTERED_ENCOMPASSING_CHIMERA_BEDPE_FILE = "encompassing_chimeras.filtered.bedpe"

EXON_JUNCTION_TRIM_BP = 10
JUNC_REF_FASTA_FILE = "spanning_juncs.fa"
JUNC_REF_MAP_FILE = "spanning_juncs.txt"
JUNC_BOWTIE_INDEX = "spanning_juncs"
JUNC_BOWTIE_INDEX_FILE = "spanning_juncs.1.ebwt"
JUNC_READS_BAM_FILE = "aligned_junc_reads.bam"
RAW_CHIMERA_BEDPE_FILE = "chimeras.raw.bedpe"
CHIMERA_BEDPE_FILE = "chimeras.bedpe"
RANKED_CHIMERA_BEDPE_FILE = "chimeras.ranked.bedpe"

# in-place XML prettyprint formatter
def indent(elem, level=0):
    i = "\n" + level*"  "
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = i + "  "
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
        for elem in elem:
            indent(elem, level+1)
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
    else:
        if level and (not elem.tail or not elem.tail.strip()):
            elem.tail = i