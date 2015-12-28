# Introduction #

The `chimerascan` pipeline is a software tool to search high-throughput RNA sequencing (RNA-seq) data for transcription of chimeric genes.  We define a _chimera_ as a single transcript that is produced from two independent genes.  Chimeras can occur due to aberrant or "noisy" splicing between nearby genes on the same chromosome.  These types of chimeras are known as _read-throughs_.  Chimeras between distant genes, or genes on different chromosomes, could possibly occur due to other splicing-related mechanisms.  In our experience, these types of _trans-splicing_ events are rare.

Chimeras can also arise due to mutations and chromosomal aberrantions in genomic DNA. Chromosomal deletions, amplifications, and rearrangements can bring distant genes close together.  When these genes are transcribed, they may splice together and form a chimera.  The term _gene fusion_ can be used to describe chimeras that occur due to genomic aberrations.

# Software requirements #

`Chimerascan` was written with minimal software dependencies and additional requirements.  We hope that users will find it easy to install and use.

## Linux/Unix ##

`chimerascan` makes some system calls that assume a standard linux or unix environment.  Specifically, it calls the 'sort' command to sort certain intermediate files.  Your linux/mac distribution must have the 'sort' command in the system path.  We tested `chimerascan` on a 64-bit Redhat Enterprise Linux 5 platform.

## python 2.6+ ##

The `chimerascan` program is written in python 2.6+.  If you do not have python2.6+ installed on your system, visit http://python.org and install it.

## The Jinja2 python package (optional) ##

The Jinja2 package is used to produce HTML output files.  This feature is optional but recommended for easy readability. You can easily install this package from the [Jinja2](http://jinja.pocoo.org/docs/intro/#installation) website.

## bowtie ##

`chimerascan` uses the **bowtie** alignment tool as the kernel for much of its searching process.  You can download and install bowtie from its Sourceforge page at (http://bowtie-bio.sourceforge.net/).  We tested `chimerascan` using bowtie 0.12.7, and cannot confirm its compatibility with older (or newer) versions.  However, `chimerascan` is likely to work with any of bowtie that supports the SAM output format.  It would be preferable to place the bowtie program into your system path, so that merely typing `bowtie` at the linux command prompt suffices to run the tool.

# Installing #

Download `chimerascan` from the Downloads page of this site.

Unpack the `.tar.gz` file as follows:
```
$ tar -zxvf chimerascan-<version>.tar.gz
```

Next, change to the `chimerascan-<version>` directory that was just created.  Run the following python command to build the tool.  This will compile any C extension modules.

```
$ python setup.py build
```

Then run the install command as follows:

```
$ python setup.py install [--prefix <INSTALL_DIR]
```

If you specified `--prefix`, `chimerascan` will be installed to a custom directory. The scripts reside in `INSTALL_DIR/bin` and the library modules reside in `INSTALL_DIR/lib/<python-version>/site-packages/chimerascan`.

If you installed `chimerascan` to a custom location, you will need to set your PYTHONPATH variable to point to the installation directory. In the `bash` shell, this can be done as follows:

```
$ export PYTHONPATH=<INSTALL_DIR>/lib/<python-version>/site-packages:$PYTHONPATH
```

You can add this to one of your login scripts (e.g. `.bashrc`) to automatically add chimerascan libraries to your `PYTHONPATH` every time you open a new terminal.

To test that the chimerascan libraries are in your `PYTHONPATH` execute the following from a terminal prompt:
```
$ python
Python 2.6.X 
Type "help", "copyright", "credits" or "license" for more information.

>>> import chimerascan
>>> chimerascan.__version__
'0.4.0'
>>> exit()
```

If this finishes without errors, you should be ready to run the `chimerascan` indexer.

# Building an index for alignment using the **chimerascan** indexer #

`chimerascan` works by aligning short sequence reads to genome and transcriptome reference sequences.  We have tested `chimerascan` extensively on the human genome using version (hg19) and the UCSC known genes, but it should also work with other genomes and gene annotation databases.

To build an index for chimerascan, you will need **both** a **genome** and a **transcriptome** reference.  The **genome** reference must be in FASTA format, and the **transcriptome** reference must be prepared in a custom format similar to BED.

## Downloading a genome reference from UCSC ##

You can easily download a genome reference from UCSC using either the `rsync` or `wget` commands.  For example, to download the latest version of the human genome (hg19), use the following command:
```
$ rsync -avzP rsync://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz .
```

Visit the UCSC download site at [ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips](ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips) for more instructions (see the README.txt) file.

## Uncompressing and preparing the genome reference ##

You must uncompress and concatenate the individual chromosome FASTA files into a single genome FASTA file before running the `chimerascan` indexer.

To uncompress, type:
```
$ tar -zxvf chromFa.tar.gz
```

We usually want to use the canonical chromosome references and ignore the haplotypes and unassembled FASTA files.  For the human genome (hg19), the following command can be used:

```
$ cat chr?.fa chr??.fa > hg19.fa
```

This will concatenate all chromosomes with one or two characters following the 'chr' file prefix.

### Obtaining transcriptome annotations from this site ###

A version of the UCSC known genes for human genome version hg19 has been posted in the Downloads section for convenience.  If you are processing human data and would be satisfied with using the UCSC genes, you can download and unzip this file using `gunzip` and proceed to the next section.

### Downloading transcriptome annotations from UCSC ###

  * Visit the UCSC Tables page at http://genome.ucsc.edu/cgi-bin/hgTables?command=start
  * Set the **clade**, **genome**, and **assembly** flags to match the genome you downloaded
  * In the **group** list select "Genes and Gene Prediction Tracks"
  * From the **track** list select "UCSC Genes" (or another list of your choice)
  * From **output format** select "selected fields from primary and related tables"
  * In the **file type returned** box enter a file name of your choice

**Click 'get output'**.  You will be redirected to a second page that allows you to select tables fields to include in the output.

In the upper table (mine reads "Select Fields from hg19.knownGene") select all fields **EXCEPT** proteinID and alignID.

In the lower table (mine reads "hg19.kgXref fields") select "geneSymbol" and no other fields.

**Click 'get output'** again.  Save the file and note the path.

If you downloaded a compressed (gzipped) transcriptome annotation file, uncompress it as follows:

```
$ gunzip <myfile.gz>
```

### Building a custom transcriptome reference ###

The transcriptome annotation reference follows a format similar to the BED format in UCSC.  Support for BED, GFF, and other formats is planned for a future release.  In the meantime, you will need to set up your gene annotations in the following format:

| **Gene ID** | **chrom** | **strand** | **txStart** | **txEnd** | **cdsStart** | **cdsEnd** | **exonCount** | **exonStarts** | **exonEnds** | **geneSymbol**|
|:------------|:----------|:-----------|:------------|:----------|:-------------|:-----------|:--------------|:---------------|:-------------|:|
| uc002gii.1  | chr17     | -          | 7571719     | 7578811   | 7572926      | 7578533    | 7             | 7571719,7573926,7576852,7577018,7577498,7578176,7578370, | 7573008,7574033,7576926,7577155,7577608,7578289,7578811, | TP53 |

The file should be in tabular text format with fields separated by the tab character ('\t').

# Running the chimerascan indexer #

Now that you have downloaded and prepared the genome and transcriptome references you can now run the chimerascan indexer.  The indexer requires that you specify the paths to the genome FASTA file, the gene annotation file, and an output directory.  The indexer uses the `bowtie-build` indexing scheme from bowtie to build the index, so you must specify the path to this program as well.

```
$ python <INSTALL_DIR>/bin/chimerascan_index.py -h
Usage: chimerascan_index.py [options] <reference_genome.fa> <gene_models.txt> <index_output_dir>

Options:
  -h, --help            show this help message and exit
  --bowtie-dir=BOWTIE_DIR
                        Path to 'bowtie-build' program (by default, expects
                        the 'bowtie-build' binary to be in current PATH)

$ python <INSTALL_DIR>/bin/chimerascan_index.py hg19.fa genes.txt myindexdir
```

This command will take several hours for large genomes.  Once the index is built, you will not need to rebuild it.

# Running #

If you completed the previous steps successfully, you can now proceed to run `chimerascan`!  See the [Running](Running.md) page for more instructions.