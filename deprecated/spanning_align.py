'''
Created on Nov 5, 2010

@author: mkiyer
'''
import logging
import argparse
import subprocess
import tempfile
import os

def get_read_length(fastq_file):
    f = open(fastq_file)
    f.next()
    seq = f.next().strip()
    f.close()
    return len(seq)

def setup_bowtie(output_sam_file, fastq_file, fastq_format, seed_length, 
                 mismatches, num_threads, bowtie_bin, bowtie_index):
    # get the read length to determine how much trimming is needed
    args = [bowtie_bin, "-q",
            "-p", str(num_threads),
            "--tryhard",
            "--%s" % fastq_format,
            "-l", str(seed_length),
            "-n", str(mismatches),
            "-a"]
    args += [bowtie_index, fastq_file, output_sam_file]
    return args

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")    
    parser = argparse.ArgumentParser()
    parser.add_argument("--bowtie-bin", dest="bowtie_bin")
    parser.add_argument("--bowtie-index", dest="bowtie_index")
    parser.add_argument("--processors", type=int, dest="processors", default=1)
    parser.add_argument("--seed-length", type=int, dest="seed_length", default=28)
    #parser.add_argument("--multihits", type=int, dest="multihits", default=40)
    parser.add_argument("--mismatches", type=int, dest="mismatches", default=2)
    parser.add_argument("--quals", dest="fastq_format")
    parser.add_argument("fastq_file")
    parser.add_argument("output_file")
    options = parser.parse_args()
    # setup bowtie
    args = setup_bowtie(options.output_file, options.fastq_file, options.fastq_format,
                        options.seed_length, options.mismatches, options.processors, 
                        options.bowtie_bin, options.bowtie_index)
    logging.debug("Running Bowtie with args: %s" % args)
    subprocess.call(args)

if __name__ == '__main__': main()



#def main():
#    logging.basicConfig(level=logging.DEBUG,
#                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")    
#    parser = argparse.ArgumentParser()
#    parser.add_argument("--tophat-bin", dest="tophat_bin")
#    parser.add_argument("--bowtie-index", dest="bowtie_index")
#    parser.add_argument("--processors", type=int, dest="processors", default=1)
#    parser.add_argument("--multihits", type=int, dest="multihits", default=40)
#    parser.add_argument("--mismatches", type=int, dest="mismatches", default=2)
#    parser.add_argument("--quals", dest="fastq_format")
#    parser.add_argument("fastq_file")
#    parser.add_argument("output_dir")
#    options = parser.parse_args()
#    # setup tophat
#    tophat_kwargs = {"--num-threads": options.processors,
#                     "--max-multihits": options.multihits,
#                     "--segment-mismatches": options.mismatches}
#    args = setup_tophat(options.tophat_bin, options.bowtie_index, 
#                        [options.fastq_file], options.fastq_format,
#                        options.output_dir, **tophat_kwargs)
#    logging.debug("Running Tophat with args: %s" % args)
#    subprocess.call(args)
#    
#def write_bogus_juncs_file(fh, read_length):
#    print >>fh, "BLABBY0000001\t%d\t%d\t+" % (read_length - 25, read_length + 25)
#
#def setup_tophat(tophat_bin, bowtie_index, fastq_file, fastq_format,
#                 output_dir, **kwargs):
#    fd,tmp_juncs_file = tempfile.mkstemp(suffix=".juncs", prefix="tmp", dir=output_dir)
#    os.close(fd)
#    write_bogus_juncs_file(open(tmp_juncs_file, "w"), get_read_length(fastq_file))
#    args = [tophat_bin,
#            "-o", output_dir,
#            "--no-novel-juncs",
#            "--raw-juncs", tmp_juncs_file]
#    if fastq_format is not None:
#        args.append("--%s" % (fastq_format))
#    for k,v in kwargs.iteritems():
#        args.extend(["%s" % k, v])
#    args.append(bowtie_index)
#    args.append(fastq_file)
#    return map(str, args)