## Whats New ##

Chimerascan version 0.4.5 (released February 25, 2012) is primarily a bug-fix release.
  * Fixed a bug where genes spanning large genomic distances could effectively "mask" gene fusions occuring within their intronic regions.
  * Fixed a bug where some fusions were inappropriately reported as "read-through" when in fact they were the result of genomic structural alterations.  These fusions will now be reported as "Adjacent\_Complex" fusions.
  * Added support for GTF gene model files through the new tools script "gtf\_to\_genepred.py"
  * Bowtie now ignores quality scores by default and uses "-v" to align reads in all stages
  * Added Ensembl GRCH37 version 65 genes to the Downloads page

Chimerascan version 0.4.3-1 (released October 28, 2011) fixes a bug that introduced false positives whenever long reads (>40bp) were used.  Now longer reads should perform extremely well.

Chimerascan version 0.4.3 (released Aug 2, 2011) includes a number of improvements and is recommended for all users
  * Support for compressed FASTQ files (gzip, bzip2, and ZIP)
  * A new module was added to detect and remove erroneous chimeras arising from homologous gene families such as HLA
  * An updated filter now excludes chimeras whenever the number of chimeric reads is much lower than the number of wild-type reads spanning the breakpoint
  * A new module was added to resolve ambiguously mapping reads based on information about the chimera coverage, insert size distribution, and mapping quality
  * Engineering improvements have reduced memory requirements at several stages of the analysis
  * The output format has been converted back to BEDPE as per user requests
  * An improved HTML report is now produced

### Introduction ###

> Recurrent gene fusions (a.k.a. chimeras) are a prevalent class of mutations that can produce functional transcripts that contribute to cancer progression.  Recent advanced in high-throughput sequencing technologies have enabled reliable gene fusion discovery (1,2).  **chimerascan** is a software package that detects gene fusions in paired-end RNA sequencing (RNA-Seq) datasets.

#### References ####
1. Maher, C.A., et al. Transcriptome sequencing to detect gene fusions in cancer. Nature 458, 97-101 (2009).

2. Maher, C.A., et al. Chimeric transcript discovery by paired-end transcriptome sequencing. Proceedings of the National Academy of Sciences of the United States of America 106, 12353-12358 (2009).