# Before getting started #

  1. Ensure that your chimerascan installation is complete (see [Installation](Installation.md)).
  1. Build an index using the chimerascan indexer (see [Installation](Installation.md)).

# Know your data #

  1. **NEW in version 0.4.3+** Your FASTQ files can be in plain text or compressed.  `chimerascan` interprets FASTQ files that end in ".gz" as `gzip` files, ".bz2" as `bzip2` files, and ".zip" as `zip` files.  Any other file extension will be treated as plain text.
  1. The quality scores of illumina FASTQ files can vary.  Check your pipeline version to learn what quality score format your data has. Specify `sanger`, `solexa`, or `illumina` for the `--quals` parameter to choose the appropriate mode when running `chimerascan`.

# Your first run #

## Finding the TMPRSS2-ERG gene fusion in prostate cancer ##

We have tested `chimerascan` using the UCSC genome version hg19 and the corresponding UCSC known gene transcriptome annotations.  For the following example, we will assume your created an index (using the `chimerascan_index.py` script) with an output directory `<INDEX_DIR>`.  `Chimerascan` comes packages with some handy test FASTQ files that we prepared from known chimeras.  These FASTQ files contain _solexa_ type quality scores. These are located in the `/tests` subdirectory within the chimerascan download directory. You should first try chimerascan on some of these to ensure that your setup is working.  Note that the test reads provided are from **human** and only compatible with indexes build from the human genome.

  1. Assume you downloaded chimerascan to the directory `DOWNLOAD_DIR`
  1. Assume you built the chimerascan index in the directory `INDEX_DIR`

To test chimerascan on the TMPRSS2-ERG discordant reads with the default parameters, run the following command:

```
python <INSTALL_DIR>/bin/chimerascan_run.py -v --quals solexa <INDEX_DIR> 
  <DOWNLOAD_DIR>/tests/vcap_pe_53bp/TMPRSS2-ERG_1.fq 
  <DOWNLOAD_DIR>/tests/vcap_pe_53bp/TMPRSS2-ERG_2.fq 
  t2erg_output_dir
```

For example, if I downloaded chimerascan to `/home/john/downloads/chimerascan`, and installed chimerascan to `/home/john/chimerascan` and built and index at '/home/john/chimerascan\_hg19\_ucsc\_index', I could run chimerascan as follows:

```
$ /home/john/chimerascan/bin/chimerascan_run.py /home/john/chimerascan_hg19_ucsc_index /home/john/downloads/chimerascan/tests/vcap_pe_53bp/TMPRSS2-ERG_1.fq /home/john/downloads/chimerascan/tests/vcap_pe_53bp/TMPRSS2-ERG_2.fq t2erg_output_dir
```

## Logging messages ##

Chimerascan produces periodic logging messages during runtime.  Adding the `--verbose` option on the command line will show additional debugging messages. If all goes well, you will see the following message upon completion:

`2011-02-09 22:04:09,907 - root - INFO - Finished run.`

A successful chimerascan run will produce several output files.  Currently, the key output file is a tabular text file named **`chimeras.bedpe`**.

You can now proceed to look at the output.

## Output ##

The `chimeras.bedpe` file contains information about the chromosomal regions, transcript ids, genes, and statistics for each chimera.  The file adapts to the BEDPE format for representing paired-intervals (courtesy Aaron Quinlan and the BEDTools project). A full description of the output format is provided below:

| **Column#** | **Name** | **Example** | **Description** |
|:------------|:---------|:------------|:----------------|
| 1           | 5' chromosome | chr21       | Chromosome of 5' partner |
| 2           | 5' start | 42870000    | Genomic start of 5' partner |
| 3           | 5' end   | 42880007    | Genomic end of 5' partner |
| 4           | 3' chromosome | chr21       | Chromosome of 3' partner |
| 5           | 3' start | 39810000    | Genomic start of 3' partner |
| 6           | 3' end   | 39817543    | Genomic end of 3' partner |
| 7           | chimera cluster id | CLUSTER001  | Identifies a group of overlapping chimeras |
| 8           | score    | 186         | Total fragments supporting chimera |
| 9           | 5' strand | -           | Genomic strand of 5' partner |
| 10          | 3' strand | -           | Genomic strand of 3' partner |
| 11          | 5' transcript ids | uc002yzj.2:0-78 | Comma-delimited list of 5' transcript ids with coordinates in chimera |
| 12          | 3' transcript ids | uc010gnw.2:200-3000,uc010gnx.2:210-2250 | Comma-delimited list of 3' transcript ids with coordinates of chimera |
| 13          | 5' genes | TMPRSS2     | Comma-delimited list of 5' genes |
| 14          | 3' genes | ERG         | Comma-delimited list of 3' genes |
| 15          | chimera type | Intrachromosomal | describes chimera as one of several types, including "Read\_Through", "Intrachromosomal", "Interchromosomal", and others |
| 16          | distance to 3' partner | -2802774    | distance from 5' to 3' chimera, or 'None' for interchromosomal chimeras |
| 17          | total fragments | 186.0       | total fragments supporting chimera (same as column 8) |
| 18          | spanning fragments | 68          | number of fragments spanning breakpoint junction |
| 19          | unique alignment positions | 132         | number unique read alignment positions |
| 20          | 5' isoform fraction | 0.97        | chimeric fraction of total 5' transcript reads that overlap the breakpoint |
| 21          | 3' isoform fraction | 0.53        | chimeric fraction of total 3' transcript reads that overlap the breakpoint |
| 22          | breakpoint spanning reads | >4668671/1;pos=29;strand=+,GGCGGAGGCGGAGGGCGAGGGGCGGGGAGCGCCGCCTGGAGCGCGGCAGGAAG,>1378815/1;pos=30;strand=+,GCGGAGGCGGAGGGCGAGGGGCGGGGAGCGCCGCCTGGAGCGCGGCAGGAAGC | pseudo-FASTA format containing reads that span the breakpoint for this chimera |
| 23          | chimera IDs | C2846325,C2846326,C2846327 | (for debugging purposes) unique IDs of specific chimera isoforms contributing to this cluster |

## Creating an HTML page for convenient output visualization ##

A script, "chimerascan\_html\_table.py", is provided in the installation and produces an HTML web page with a sortable table of the chimeras in the output. **NOTE: This script requires the Jinja2 python package to be installed, visit http://jinja.pocoo.org/docs/intro/#installation to install it**.  The script suppresses read-through chimeras by default, unless the `--read-throughs` option is supplied on the command line. The script is run as follows:

```
$ python <INSTALL_DIR/bin/chimerascan_html_table.py -h
Usage: chimerascan_html_table.py [options] <chimeras.bedpe>

Options:
  -h, --help       show this help message and exit
  -o OUTPUT_FILE   output file [default=stdout]
  --read-throughs  include read-through chimeras in output [default=False]
```

For example:
```
python /home/john/chimerascan/bin/chimerascan_html_table.py --read-throughs -o chimeras.html ./t2erg_output_dir/chimeras.bedpe
```

The output file, **chimeras.html** will be written.  You can open it in your favorite web browser.  Javascript should be enabled to facilitate the sorting feature.

# `Chimerascan` command-line options #

**These may change in subsequent release of `chimerascan` and this wiki page may lag the software**. Please run `chimerascan_run.py` with the "-h" option to view the most recent list of options.

Many of chimerascan's options are directly passed to `Bowtie`.  Therefore, users should refer to Bowtie's manual pages for a full explanation of these parameters.

```
$ <INSTALL_DIR>/bin/chimerascan_run.py --help
Usage: chimerascan_run.py [options] [--config <config_file>  | <index> <mate1.fq> <mate2.fq> <output_dir>]

Options:
  --version             show program's version number and exit
  -h, --help            show this help message and exit
  --config=CONFIG_FILE  Path to configuration XML file
  -v, --verbose         enable verbose logging output [default=False]
  -p NUM_PROCESSORS, --processors=NUM_PROCESSORS
                        Number of processor cores to allocate to chimerascan
                        [default=2]
  --keep-tmp            DO NOT delete intermediate files after run
                        [default=True]
  --rm-tmp              Delete intermediate files after run [default=True]
  --quals=FMT           FASTQ quality score format [default=sanger]
  --library-type=LIBRARY_TYPE
                        Library type [default=fr-unstranded]
  --isize-mean=N        Mean insert size to sample from when insert size
                        distribution cannot be determined empirically
                        [default=200]
  --isize-stdev=N       Insert size standard deviation to sample from when
                        insert size distribution cannot be determined
                        empirically [default=40]

  Bowtie options:
    Adjust these options to change bowtie alignment settings

    --bowtie-path=BOWTIE_PATH
                        Path to directory containing 'bowtie' [default:
                        assumes bowtie is in the current PATH]
    --trim5=N           Trim N bases from 5' end of read
    --trim3=N           Trim N bases from 3' end of read
    --min-fragment-length=MIN_FRAGMENT_LENGTH
                        Smallest expected fragment length [default=0]
    --max-fragment-length=MAX_FRAGMENT_LENGTH
                        Largest expected fragment length (reads less than this
                        fragment length are assumed to be unspliced and
                        contiguous) [default=1000]
    --multihits=N       Ignore reads that map to more than N locations
                        [default=100]
    --initial-mismatches=N
                        Aligned reads must have <= N mismatches [default=2]
    --initial-bowtie-args=BOWTIE_ARGS
                        Additional arguments to pass to bowtie aligner (see
                        Bowtie manual for options) [default='--best --strata']
    --discord-bowtie-args=DISCORD_BOWTIE_ARGS
                        Additional arguments to pass to bowtie aligner for
                        discordant alignment phase (see Bowtie manual for
                        options) [default='--best']
    --discord-mismatches=N
                        Discordant aligned reads must have <= N mismatches
                        [default=3]
    --segment-length=N  Size of read segments during discordant alignment
                        phase [default=25]

  Filtering options:
    Adjust these options to change filtering behavior

    --homology-mismatches=N
                        Number of mismatches to tolerate at breakpoints when
                        computing homology [default=2]
    --anchor-min=N      Minimum breakpoint overlap (bp) required to call
                        spanning reads [default=4]
    --anchor-length=N   Size (bp) of anchor region where mismatch checks are
                        enforced [default=8]
    --anchor-mismatches=N
                        Number of mismatches allowed within anchor region
                        [default=0]
    --filter-unique-frags=N
                        Filter chimeras with less than N unique aligned
                        fragments [default=2]
    --filter-isize-prob=X
                        Filter chimeras when probability of observing the
                        putative insert size is less than X (0.0-1.0)
                        [default=0.01]
    --filter-isoform-fraction=X
                        Filter chimeras with expression ratio less than X
                        (0.0-1.0) relative to the wild-type transcript levels
                        [default=0.01]
    --filter-false-pos=FILTER_FALSE_POS_FILE
                        File containing known false positive chimeric
                        transcript pairs to filter out
```

## `Chimerascan` configuration file support ##

Chimerascan supports an XML configuration file containing a run configuration and all command-line options.  All of the XML fields are optional and the command line parameters will always override the XML fields.  At the beginning of each run, a file `runconfig.xml` will be written to the output directory to reflect the exact conglomeration of parameters that were used in that run.

Below is an example working configuration file that may be used as a template:

```
<chimerascan_config>
  <output_dir>vcap042test</output_dir>
  <fastq_files>
    <file mate="0">/home/john/chimerascan/tests/vcap_pe_53bp/TMPRSS2-ERG_1.fq</file>
    <file mate="1">/home/john/chimerascan/tests/vcap_pe_53bp/TMPRSS2-ERG_2.fq</file>
  </fastq_files>
  <index>/exds/projects/chimerascan/hg19_ucsc_illumina</index>
  <num_processors>4</num_processors>
  <quals>solexa</quals>
  <keep_tmp>True</keep_tmp>
  <bowtie_path />
  <bowtie_args>--best --strata</bowtie_args>
  <discord_bowtie_args>--best</discord_bowtie_args>
  <multihits>100</multihits>
  <mismatches>2</mismatches>
  <discord_mismatches>3</discord_mismatches>
  <segment_length>25</segment_length>
  <trim5>0</trim5>
  <trim3>0</trim3>
  <min_fragment_length>0</min_fragment_length>
  <max_fragment_length>1000</max_fragment_length>
  <library_type>fr-unstranded</library_type>
  <isize_mean>200</isize_mean>
  <isize_stdev>40</isize_stdev>
  <homology_mismatches>2</homology_mismatches>
  <anchor_min>4</anchor_min>
  <anchor_length>8</anchor_length>
  <anchor_mismatches>0</anchor_mismatches>
  <filter_unique_frags>2.0</filter_unique_frags>
  <filter_isize_prob>0.01</filter_isize_prob>
  <filter_isoform_fraction>0.01</filter_isoform_fraction>
  <filter_false_pos_file>hg19_bodymap_false_positive_chimeras.txt</filter_false_pos_file>
</chimerascan_config>

```