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

BASE_PROCESSORS = 2
MIN_SEGMENT_LENGTH = 20
ALIGN_INDEX = 'align_index'
BOWTIE_INDEX_FILE = 'align_index.1.ebwt'
GENE_REF_PREFIX = 'gene_'
GENE_FEATURE_FILE = "gene_features.txt"

RUNCONFIG_XML_FILE = "runconfig.xml"
ALIGNED_READS_BAM_FILE = "aligned_reads.bam"
UNALIGNED_FASTQ_PARAM = "unaligned.fq"
UNALIGNED_FASTQ_FILES = ("unaligned_1.fq", "unaligned_2.fq")
MAXMULTIMAP_FASTQ_PARAM = "maxmulti.fq"
MAXMULTIMAP_FASTQ_FILES = ("maxmulti_1.fq", "maxmulti_2.fq")
ISIZE_DIST_FILE = "isize_dist.txt"
ISIZE_MAX_SAMPLES = 1e6
ISIZE_NUM_STDEVS = 3
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


