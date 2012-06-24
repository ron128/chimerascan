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
# return codes
JOB_SUCCESS = 0
JOB_ERROR = 1

# binary program names
BOWTIE2_BIN = "bowtie2"
BOWTIE2_BUILD_BIN = "bowtie2-build"
BOWTIE2_INSPECT_BIN = "bowtie2-inspect"

# constants for index
TRANSCRIPTOME_INDEX = 'transcriptome'
TRANSCRIPTOME_FASTA_FILE = 'transcriptome.fa'
TRANSCRIPT_FEATURE_FILE = 'transcripts.txt'
GENOME_INDEX = 'genome'
GENOME_FASTA_FILE = 'genome.fa'
BOWTIE2_INDEX_FILE_EXTS = ('.1.bt2', '.2.bt2', '.3.bt2', '.4.bt2', '.rev.1.bt2', '.rev.2.bt2')
GENOME_BOWTIE2_FILES = ((GENOME_INDEX + x) for x in BOWTIE2_INDEX_FILE_EXTS)  
TRANSCRIPTOME_BOWTIE2_FILES = ((TRANSCRIPTOME_INDEX + x) for x in BOWTIE2_INDEX_FILE_EXTS) 
MAX_MULTIMAPPING_FILE = 'max_multihits.txt'

# chimerascan subdirectories
LOG_DIR = "log"
TMP_DIR = "tmp"

# defaults and constraints for run configuration
BASE_PROCESSORS = 2
MIN_SEGMENT_LENGTH = 25
DEFAULT_MIN_FRAG_LENGTH = 0
DEFAULT_MAX_FRAG_LENGTH = 1000
DEFAULT_MAX_MULTIHITS = 1
RUNCONFIG_XML_FILE = "runconfig.xml"

# output after read inspection, name conversion, and 
# quality score conversion
CONVERTED_FASTQ_PREFIX = "reads"
CONVERTED_FASTQ_FILES = tuple(CONVERTED_FASTQ_PREFIX + "_%d.fq" % (x+1) 
                              for x in xrange(2))
READ_NAME_DBM_FILE = CONVERTED_FASTQ_PREFIX + ".dbm"

# output from initial alignment
TRANSCRIPTOME_BAM_FILE = "transcriptome_reads.bam"
TRANSCRIPTOME_UNALIGNED_PATH = "transcriptome_unaligned.fq"
TRANSCRIPTOME_UNALIGNED_FASTQ_FILES = ("transcriptome_unaligned.1.fq", 
                                       "transcriptome_unaligned.2.fq")
TRANSCRIPTOME_LOG_FILE = "transcriptome_alignment.log"
SORTED_TRANSCRIPTOME_BAM_FILE = "transcriptome_reads.srt.bam"

GENOME_BAM_FILE = "genome_reads.bam"
GENOME_UNALIGNED_PATH = "genome_unaligned.fq"
GENOME_UNALIGNED_FASTQ_FILES = ("genome_unaligned.1.fq", 
                                "genome_unaligned.2.fq")
GENOME_LOG_FILE = "genome_alignment.log"

# insert size estimation parameters
ISIZE_MIN_SAMPLES = 100
ISIZE_MAX_SAMPLES = 5e6
ISIZE_DIST_FILE = "isize_dist.txt"

# interleaved segmented paired-end reads file
INTERLEAVED_FASTQ_FILE = "interleaved_reads.fq"
INTERLEAVED_TRIMMED_FASTQ_FILE = "interleaved_trimmed_reads.fq"

# output from realignment of trimmed reads
REALIGNED_BAM_FILE = "realigned_reads.bam"
REALIGNED_LOG_FILE = "realigned_reads.log"

# output for different classes of discordant reads
PAIRED_BAM_FILE = "realigned_paired_reads.bam"
DISCORDANT_BAM_FILE = "realigned_discordant_pairs.bam"
UNPAIRED_BAM_FILE = "realigned_unpaired_reads.bam"
UNMAPPED_BAM_FILE = "realigned_unmapped_reads.bam"
MULTIMAP_BAM_FILE = "realigned_multimap_reads.bam"
UNRESOLVED_BAM_FILE = "realigned_unresolved_reads.bam"

# discordant pairs BAM file
DISCORDANT_GENOME_SAM_FILE = "realigned_discordant_pairs.genome.sam"
DISCORDANT_GENOME_BAM_FILE = "realigned_discordant_pairs.genome.bam"
SORTED_DISCORDANT_GENOME_BAM_FILE = "realigned_discordant_pairs.genome.srt.bam"

# discordant clusters
DISCORDANT_CLUSTER_FILE = "discordant_clusters.txt"
DISCORDANT_CLUSTER_SHELVE_FILE = "discordant_clusters.shelve"
DISCORDANT_CLUSTER_PAIR_FILE = "discordant_cluster_pairs.txt"
SORTED_DISCORDANT_GENOME_CLUSTER_BAM_FILE = "realigned_discordant_pairs.genome.clustered.srt.bam"

# output files
UNFILTERED_OUTPUT_FILE = "chimeras.unfiltered.bedpe"
OUTPUT_FILE = "chimeras.bedpe"

## discordant reads BEDPE file
#DISCORDANT_BEDPE_FILE = "discordant_reads.bedpe"
#SORTED_DISCORDANT_BEDPE_FILE = "discordant_reads.srt.bedpe"
#
## chimera candidates with encompassing read support
#ENCOMPASSING_CHIMERA_FILE = "encompassing_chimeras.txt"
#FILTERED_ENCOMPASSING_CHIMERA_FILE = "encompassing_chimeras.filtered.txt"
#
## amount of trimming to use to stop reads from overlapping 
## exon boundaries and going into intronic space
#EXON_JUNCTION_TRIM_BP = 10
#
## number of homology mismatches in breakpoint sequences 
## to tolerate when computing homology distance
#BREAKPOINT_HOMOLOGY_MISMATCHES = 2
#BREAKPOINT_CHIMERA_FILE = "encompassing_chimeras.breakpoint_sorted.txt"
#BREAKPOINT_MAP_FILE = "breakpoints.txt"
#BREAKPOINT_FASTA_FILE = "breakpoints.fa"
#BREAKPOINT_INDEX = "breakpoints"
#BREAKPOINT_BOWTIE2_FILES = ((BREAKPOINT_INDEX + x) for x in BOWTIE2_INDEX_FILE_EXTS) 
#
## reads to remap to breakpoint junction index
#ENCOMP_SPANNING_FASTQ_FILE = "encomp_spanning_reads.fq"
#SINGLE_MAPPED_BAM_FILE = "singlemap_reads.srt.bam"
#SINGLEMAP_SPANNING_FASTQ_FILE = "singlemap_spanning_reads.fq"
#UNALIGNED_SPANNING_FASTQ_FILE = "unaligned_spanning_reads.fq"
#
## results of aligning reads to breakpoint index
#ENCOMP_SPANNING_BAM_FILE = "encomp_spanning_reads.bam"
#SORTED_ENCOMP_SPANNING_BAM_FILE = "encomp_spanning_reads.srt.bam"
#SINGLEMAP_SPANNING_BAM_FILE = "singlemap_spanning_reads.bam"
#SORTED_SINGLEMAP_SPANNING_BAM_FILE = "singlemap_spanning_reads.srt.bam"
#UNALIGNED_SPANNING_BAM_FILE = "unaligned_spanning_reads.bam"
#SORTED_UNALIGNED_SPANNING_BAM_FILE = "unaligned_spanning_reads.srt.bam"
#
## results of merging spanning reads into chimera nominations
#SPANNING_CHIMERA_FILE = "spanning_chimeras.txt"
## results of resolving ambiguous reads
#RESOLVED_SPANNING_CHIMERA_FILE = "spanning_chimeras.resolved.txt"
## results of filtering chimeras
#FILTERED_CHIMERA_FILE = "spanning_chimeras.resolved.filtered.txt"
#HOMOLOG_FILTERED_CHIMERA_FILE = "spanning_chimeras.resolved.filtered.homolog.txt"
#BEST_FILTERED_CHIMERA_FILE = "spanning_chimeras.resolved.filtered.homolog.best_isoform.txt"
## output file
#CHIMERA_OUTPUT_FILE = "chimeras.bedpe"