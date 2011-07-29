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

# constants for index
ALIGN_INDEX = 'align_index'
ALIGN_INDEX_FASTA_FILE = 'align_index.fa'
BOWTIE_INDEX_FILE = 'align_index.1.ebwt'
FRAG_SIZE_INDEX = "frag_size_index"
FRAG_SIZE_INDEX_FILE = "frag_size_index.1.ebwt"
GENE_REF_PREFIX = 'gene_'
GENE_FEATURE_FILE = "gene_features.txt"
RAW_JUNCS_FILE = "known_juncs.txt"

# chimerascan subdirectories
LOG_DIR = "log"
TMP_DIR = "tmp"

# constraints for run configuration
BASE_PROCESSORS = 2
MIN_SEGMENT_LENGTH = 20
RUNCONFIG_XML_FILE = "runconfig.xml"

# output after read inspection, name conversion, and 
# quality score conversion
CONVERTED_FASTQ_PREFIX = "reads"
CONVERTED_FASTQ_FILES = tuple(CONVERTED_FASTQ_PREFIX + "_%d.fq" % (x+1) 
                              for x in xrange(2))

# output from initial bowtie alignment
ALIGNED_READS_BAM_FILE = "aligned_reads.bam"
UNALIGNED_FASTQ_PARAM = "unaligned.fq"
UNALIGNED_FASTQ_FILES = ("unaligned_1.fq", "unaligned_2.fq")
MAXMULTIMAP_FASTQ_PARAM = "maxmulti.fq"
MAXMULTIMAP_FASTQ_FILES = ("maxmulti_1.fq", "maxmulti_2.fq")

# sorted aligned reads bam file
SORTED_ALIGNED_READS_BAM_FILE = "sorted_aligned_reads.bam"

# insert size estimation parameters
ISIZE_MIN_SAMPLES = 100
ISIZE_MAX_SAMPLES = 1e6
ISIZE_DIST_FILE = "isize_dist.txt"

# output from realignment of trimmed reads
REALIGNED_BAM_FILE = "realigned_reads.bam"

# output for different classes of discordant reads
GENE_PAIRED_BAM_FILE = "gene_paired_reads.bam"
GENOME_PAIRED_BAM_FILE = "genome_paired_reads.bam"
REALIGNED_UNMAPPED_BAM_FILE = "unmapped_reads.bam"
REALIGNED_COMPLEX_BAM_FILE = "complex_reads.bam"

# discordant reads BEDPE file
DISCORDANT_BEDPE_FILE = "discordant_reads.bedpe"
SORTED_DISCORDANT_BEDPE_FILE = "discordant_reads.srt.bedpe"

# chimera candidates with encompassing read support
ENCOMPASSING_CHIMERA_FILE = "encompassing_chimeras.txt"

# amount of trimming to use to stop reads from overlapping 
# exon boundaries and going into intronic space
EXON_JUNCTION_TRIM_BP = 10

# number of homology mismatches in breakpoint sequences 
# to tolerate when computing homology distance
BREAKPOINT_HOMOLOGY_MISMATCHES = 2
BREAKPOINT_CHIMERA_FILE = "encompassing_chimeras.breakpoint_sorted.txt"
BREAKPOINT_MAP_FILE = "breakpoints.txt"
BREAKPOINT_FASTA_FILE = "breakpoints.fa"
BREAKPOINT_BOWTIE_INDEX = "breakpoints"
BREAKPOINT_BOWTIE_INDEX_FILE = "breakpoints.1.ebwt"

# reads to remap to breakpoint junction index
ENCOMP_SPANNING_FASTQ_FILE = "encomp_spanning_reads.fq"
SINGLEMAP_SPANNING_FASTQ_FILE = "singlemap_spanning_reads.fq"
UNALIGNED_SPANNING_FASTQ_FILE = "unaligned_spanning_reads.fq"
# results of aligning reads to breakpoint index
ENCOMP_SPANNING_BAM_FILE = "encomp_spanning_reads.bam"
SORTED_ENCOMP_SPANNING_BAM_FILE = "encomp_spanning_reads.srt.bam"
SINGLEMAP_SPANNING_BAM_FILE = "singlemap_spanning_reads.bam"
SORTED_SINGLEMAP_SPANNING_BAM_FILE = "singlemap_spanning_reads.srt.bam"
UNALIGNED_SPANNING_BAM_FILE = "unaligned_spanning_reads.bam"
SORTED_UNALIGNED_SPANNING_BAM_FILE = "unaligned_spanning_reads.srt.bam"

# results of merging spanning reads into chimera nominations
SPANNING_CHIMERA_FILE = "spanning_chimeras.txt"
# results of resolving ambiguous reads
RESOLVED_SPANNING_CHIMERA_FILE = "resolved_spanning_chimeras.txt"
# results of filtering chimeras
FILTERED_CHIMERA_FILE = "spanning_chimeras.filtered.txt"
# output file
CHIMERA_OUTPUT_FILE = "chimeras.bedpe"

