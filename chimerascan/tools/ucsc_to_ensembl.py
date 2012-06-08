'''
Created on Nov 1, 2011

@author: mkiyer
'''
import argparse
import sys

chrom_names = [["chr1", "1"],
               ["chr2", "2"],
               ["chr3", "3"],
               ["chr4", "4"],
               ["chr5", "5"],
               ["chr6", "6"],
               ["chr7", "7"],
               ["chr8", "8"],
               ["chr9", "9"],
               ["chr10", "10"],
               ["chr11", "11"],
               ["chr12", "12"],
               ["chr13", "13"],
               ["chr14", "14"],
               ["chr15", "15"],
               ["chr16", "16"],
               ["chr17", "17"],
               ["chr18", "18"],
               ["chr19", "19"],
               ["chr20", "20"],
               ["chr21", "21"],
               ["chr22", "22"],
               ["chrX", "X"],
               ["chrY", "Y"],
               ["chrM", "MT"]]

ucsc_to_ensembl_dict = dict(chrom_names)
ensembl_to_ucsc_dict = dict((v,k) for k,v in ucsc_to_ensembl_dict.iteritems())

parser = argparse.ArgumentParser()
parser.add_argument("-r", action="store_true", default=False, dest="reverse")
parser.add_argument("--col", dest="colnum", type=int, default=0)
parser.add_argument("filename")
args = parser.parse_args()
colnum = max(0, args.colnum - 1)

if args.reverse:
    mydict = ensembl_to_ucsc_dict
else:
    mydict = ucsc_to_ensembl_dict

for line in open(args.filename):
    if not line:
        continue
    line = line.strip()
    if line.startswith("track") or line.startswith("#"):
        print line
    else: 
        fields = line.strip().split('\t')    
        if fields[colnum] not in mydict:
            print >>sys.stderr, "skipped line %s" % line
            continue
        fields[colnum] = mydict[fields[colnum]]
        print '\t'.join(fields)
