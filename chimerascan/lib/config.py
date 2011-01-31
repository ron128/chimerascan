'''
Created on Jan 5, 2011

@author: mkiyer
'''
'''
Created on Nov 16, 2010

@author: mkiyer
'''
JOB_SUCCESS = 0
JOB_ERROR = 1

BASE_PROCESSORS = 2
MIN_SEGMENT_LENGTH = 20
ALIGN_INDEX = 'align_index'
BOWTIE_INDEX_FILE = 'align_index.1.ebwt'
GENE_REF_PREFIX = 'gene_'
GENE_FEATURE_FILE = "gene_features.txt"

ALIGNED_READS_BAM_FILE = "aligned_reads.bam"
UNALIGNED_FASTQ_PARAM = "unaligned.fq"
UNALIGNED_FASTQ_FILES = ("unaligned_1.fq", "unaligned_2.fq")
MAXMULTIMAP_FASTQ_PARAM = "maxmulti.fq"
MAXMULTIMAP_FASTQ_FILES = ("maxmulti_1.fq", "maxmulti_2.fq")
ISIZE_STATS_FILE = "isize_stats.txt"
DISCORDANT_BAM_FILE = "discordant_reads.bam"
DISCORDANT_PAIRED_BAM_FILE = "discordant_reads_paired.bam"
DISCORDANT_BEDPE_FILE = "discordant_reads.bedpe"
EXTENDED_DISCORDANT_BEDPE_FILE = "discordant_reads.extended.bedpe"
SORTED_DISCORDANT_BEDPE_FILE = "discordant_reads.srt.bedpe"
SPANNING_FASTQ_FILE = "putative_spanning_reads.fq"
ENCOMPASSING_CHIMERA_BEDPE_FILE = "encompassing_chimeras.bedpe"
EXON_JUNCTION_TRIM_BP = 10
JUNC_REF_FASTA_FILE = "spanning_juncs.fa"
JUNC_REF_MAP_FILE = "spanning_juncs.txt"
JUNC_BOWTIE_INDEX = "spanning_juncs"
JUNC_READS_BAM_FILE = "aligned_junc_reads.bam"
RAW_CHIMERA_BEDPE_FILE = "chimeras.raw.bedpe"
CHIMERA_BEDPE_FILE = "chimeras.bedpe"

