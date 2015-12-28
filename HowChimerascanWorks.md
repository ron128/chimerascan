# How `chimerascan` works #

The following is a step-by-step description of the chimerascan workflow.

## Step 0: Create alignment index ##

The `chimerascan_index.py` script included with `chimerascan` is used to prepare an index for alignment (see [Installation](Installation.md)).

## Step 1: Prepare reads for alignment ##

`chimerascan` parses FASTQ files containing the original read sequences and performs the following: 1) converts all quality scores to Sanger format (Phred + 33) and 2) converts the `qname` for the reads from an arbitrarily long string to a number.  Renaming the reads guarantees that the paired-end naming convention that uses "/1" and "/2" to denote read 1 and read 2 is followed, and also reduces the storage requirements of intermediate steps. This process essentially creates a temporary copy of the original reads.

## Step 2: Align paired-end reads ##

Paired-end reads are aligned to a combined genome and transcriptome reference.  By design this creates many redundant alignments for each read because a single read could align to both the genome and to different transcript isoforms.  Therefore, the alignment process allows for multiple mappings per read to accommodate this redundancy (`--multimaps, default: 40`).

This initial alignment is performed in paired-end mode; both reads in a pair must align within a distance range specified by the user with the `--min-fragment-length` and `--max-fragment-length` flags. The default settings use a fragment length range `0-1000bp`.

There are a number of other command-line options that are used to control the behavior of Bowtie at this initial paired-end alignment stage.  These are:

  1. Read trimming (`--trim5` and `--trim3`) trims a specified number of bases from the 5' or 3' end of all reads.
  1. Number of mismatches (`--mismatches`) tolerated in alignments
  1. Additional arguments (`--bowtie-args`) can be directly passed to Bowtie.  **NOTE** It is not recommended to use this unless the user is absolutely sure about it!


## Step 3: Create a sorted/indexed BAM file ##

A sorted, indexed BAM file is created to enable fast lookup of the original read alignments by genomic coordinates. This employs the `pysam` package (bundled within `chimerascan`) to interact with the C function in the `samtools` library (http://code.google.com/p/pysam/).

This enables fast lookup of reads for downstream analyses such as gene expression. This process can take a considerable amount of time.

## Step 4: Estimate fragment size distribution ##

`chimerascan` parses the aligned reads and computes the empirical distribution of fragment sizes.  Only uniquely mapping reads are used to sample the fragment size distribution. The fragment size distribution is used in future steps to help localize fusion breakpoints.

## Step 5: Realign initially unmapped reads ##

All of the read pairs that successfully mapped in the initial alignment step are called **concordant reads**, where both reads in the pair localize to a specific transcript or genomic region and do not provide evidence for the existence of chimeras.

Some of read pairs that could not be mapped to either the genome or transcriptome may be **discordant** e.g. the individual reads in the pair align to different transcripts or distance regions of the genome.

In the second alignment step, all of the initially unmapped reads are treated as single reads and realigned.  Additionally, the reads are trimmed such that only the sequences at the ends of the fragment are aligned.  The size of the trimmed segment can be specified by the user using the `--segment-length` option (default=25bp). Setting the `--segment-length` to a large number will increase the specificity of the mapping process (less multi-mapping reads) but will decrease sensitivity to detect chimeras, since a larger percentage of reads will contain fusion breakpoints and fail alignment. If you suspect that your library preparation produced short DNA fragments sequenced at long read lengths, then a considerable fraction of the read pairs could be overlapping.  In this case try to keep `--segment-length` small.

A number of command-line options can be changed to control the behavior of Bowtie at the realignment stage:
  1. The number of mismatches tolerated (`--discord-mismatches`, default 3).  It is **recommended** to allow for 3 mismatches here (the maximum supported by Bowtie) as this will reduce the number of discordant reads arising from homologous genes.
  1. Additional Bowtie arguments (`--discord-bowtie-args`).  **NOTE** do not change this option unless you are sure about it.

## Step 6: Discover discordant reads ##

The realigned reads are searched for evidence of discordant reads.  A discordant read is currently defined as follows:

  1. The fragment does not align to the genome within user-specified fragment size range (default: 1000bp)
  1. The pair does not align to a single transcript
  1. The pair does not align to different transcripts that share exonic sequences on the same strand.  In otherwords, reads that map to different isoforms of the same gene are excluded, as it is not the goal of chimerascan to discover unannotated splicing patterns within a single gene.

Discordant reads are stored in an intermediate file format called BEDPE (developed by Aaron Quinlan, http://code.google.com/p/bedtools).  These reads are sorted by reference name and position.

## Step 7: Nominate chimeras ##

The fragment size distribution is used to predict an optimal fusion breakpoint location for each discordant pair along the 5' and 3' transcripts.  All reads that share the same breakpoint are aggregated together to form putative chimeras.

## Step 8: Extract chimeric breakpoint sequences ##

The upstream and downstream sequences surrounding each chimeric breakpoint are extracted from the genome FASTA file.  The `bowtie-build` indexer is used to create a new alignment index of these breakpoint sequences.

## Step 9: Nominate reads that could span breakpoints ##

Fragments that span chimera breakpoints can be divided into two classes: 1) discordant read pairs that encompass the breakpoint and aligned successfully after being trimmed and 2) read pairs that failed to align during the realignment step because one or both of the reads span a breakpoint. Reads that meet either of these criteria are written to a new FASTQ file for alignment against the breakpoint sequence index.

## Step 10: Align against breakpoint sequence database ##

The `bowtie` aligner is used to align all candidate breakpoint spanning reads to the breakpoint database from step #8.

## Step 11: Assess breakpoint spanning alignment results ##

Reads that aligned to the breakpoint sequence database are considered _spanning_ a breakpoint according to the following criteria:

  1. The alignment overlaps the breakpoint by at least `--anchor-min` bases (default: 4) **plus** the number bases of homologous sequence between the chimeric transcript and the wild-type transcript.
  1. No more than `--anchor-mismatches` mismatches are found within the first `--anchor-length` bases of the alignment.

## Step 12: Filter chimeras ##

The chimeras are passed through a number of filters in order to remove erroneous artifacts.  These include:

  1. Chimeras with very low coverage (specified with `--filter-unique-frags`, default 2) may have arisen from ligation artifacts during library preparation
  1. Chimeric transcript expression may be much lower than the expression of one or both of the wild-type transcripts in the sample (specified with `--filter-isoform-fraction`, default 0.10).
  1. Reads that support chimeras may not agree with the fragment size distribution of the library (specified with `--filter-isize-percentile`, default 99%).
  1. Users may have a list of known false positives available (specified with `--filter-false-pos` as a path to a file containing false positives)

## Step 13: Produce a text output file ##

Finally, a tab-delimited text file (currently named `chimeras.txt`) is produced in the output directory.  The fields in this file are described in the [Running](Running.md) wiki page on this site.

## (optional) Step 14: Produce an HTML web site ##

Users may run the companion script `chimerascan_html_table.py` to generate an HTML output of the chimera results for further investigation.

# Questions or concerns #

We are continuing to develop and improve chimerascan and would welcome your feedback and suggestions. Please post your issues or concerns with chimerascan on the **Issues** section of this web site.