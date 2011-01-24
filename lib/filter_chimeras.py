'''
Created on Jan 23, 2011

@author: mkiyer
'''
import collections
import logging
import operator
import subprocess
import tempfile
import os

# awk '$8 >= 2' ./mcf7_30LEJAAXX_7/chimeras.bedpe | sort -k19,19gr -k20,20gr | cut -f7-14,18-21 | less


class Chimera(object):
    @staticmethod
    def from_tabular(line):
        c = Chimera()
        fields = line.strip().split('\t')
        c.tx5p = fields[0]
        c.start5p = fields[1]
        c.end5p = fields[2]
        c.tx3p = fields[3]
        c.start3p = fields[4]
        c.end3p = fields[5]
        c.name = fields[6]
        c.gene5p, c.gene3p = fields[6].split('-', 1)
        c.encompassing_reads = int(fields[7])
        c.strand5p, c.strand3p = fields[8:10]
        c.chimera_type = fields[10]
        if fields[11] == "None":
            c.distance = None
        else:
            c.distance = int(fields[11])
        c.exons5p = fields[12]
        c.exons3p = fields[13]
        c.encompassing_ids = fields[14]
        c.encompassing_seq5p = fields[15]
        c.encompassing_seq3p = fields[16]
        c.spanning_reads = int(fields[17])
        c.encomp_and_spanning = int(fields[18])
        c.total_reads = int(fields[19])
        c.junction_hist = map(int, fields[20].split(','))
        c.spanning_seq = fields[21]
        c.spanning_ids = fields[22]
        return c
    
    def to_tabular(self):
        fields = [self.tx5p, self.start5p, self.end5p,
                  self.tx3p, self.start3p, self.end3p,
                  self.name, self.encompassing_reads,
                  self.strand5p, self.strand3p, self.chimera_type,
                  self.distance, self.exons5p, self.exons3p,
                  self.encompassing_ids, self.encompassing_seq5p, self.encompassing_seq3p,
                  self.spanning_reads, self.encomp_and_spanning, self.total_reads, 
                  ','.join(map(str, self.junction_hist)), 
                  self.spanning_seq, self.spanning_ids]
        return '\t'.join(map(str, fields))

def make_temp(base_dir, suffix=''):
    fd,name = tempfile.mkstemp(suffix=suffix, prefix='tmp', dir=base_dir)
    os.close(fd)
    return name

def filter_isoforms(input_file, tmp_dir):
    # sort by 5'/3' fusion name
    logging.debug("Sorting chimeras by 5'/3' gene name")
    tmpfile = make_temp(tmp_dir, ".bedpe")    
    fh = open(tmpfile, "w")
    subprocess.call(["sort", "-k7,7", input_file], stdout=fh)
    fh.close()
    # parse sorted file
    logging.debug("Choosing highest coverage isoforms for each gene pair")
    chimeras = []
    for line in open(tmpfile):
        c = Chimera.from_tabular(line)
        if len(chimeras) > 0 and chimeras[-1][1].name != c.name:
            # sort chimeras by total reads
            best_chimera = sorted(chimeras, key=operator.itemgetter(0))[-1][1]
            yield best_chimera
            chimeras = []
        chimeras.append((c.total_reads, c))
    if len(chimeras) > 0:
        best_chimera = sorted(chimeras, key=operator.itemgetter(0))[-1][1]
        yield best_chimera
    # remove temporary file
    os.remove(tmpfile)

def filter_overlapping_genes(chimeras):
    logging.debug("Filtering overlapping genes")
    for chimera in chimeras:
        if chimera.distance == 0:
            continue
        yield chimera

def filter_coverage(chimeras, coverage):
    logging.debug("Filtering coverage < %d" % (coverage))
    for chimera in chimeras:
        if chimera.total_reads < coverage:
            continue
        yield chimera

def filter_chimeras(input_bedpe_file, output_bedpe_file, coverage=2):
    # setup multiple filters
    it1 = filter_isoforms(input_bedpe_file, os.path.dirname(output_bedpe_file))
    it2 = filter_overlapping_genes(it1)
    it3 = filter_coverage(it2, coverage)
    # filter candidates
    logging.debug("Sorting chimeras by read coverage")
    chimeras = list(it3)
    sorted_chimeras = sorted(chimeras, key=operator.attrgetter("encomp_and_spanning", "total_reads"), reverse=True)
    outfh = open(output_bedpe_file, "w")
    for c in sorted_chimeras:
        print >>outfh, c.to_tabular()
    outfh.close()
    return

    encomp_dict = collections.defaultdict(lambda: 0)
    spanning_dict = collections.defaultdict(lambda: 0)
    total_dict = collections.defaultdict(lambda: 0)
    # profile results
    for line in open(input_bedpe_file):
        c = Chimera.from_tabular(line)
        encomp_dict[c.encompassing_reads] += 1
        spanning_dict[c.spanning_reads] += 1
        total_dict[c.total_reads] += 1    
    import numpy as np
    sorted_keys = sorted(encomp_dict)
    arr = np.array(sorted_keys[-1], dtype=np.int)
    for x in xrange(sorted_keys[-1]):
        arr[x] = encomp_dict[x]
        print x, encomp_dict[x]


def main():
    from optparse import OptionParser
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")    
    parser = OptionParser("usage: %prog [options] <chimeras.bedpe> <filtered_chimeras.bedpe>")
    parser.add_option("--overlap", action="store_true")
    parser.add_option("--encompassing-prob", type="float", dest="encomp_prob", default=0.95)
    parser.add_option("--spanning-prob", type="float", dest="spanning_prob", default=0.95)
    options, args = parser.parse_args()
    input_bedpe_file = args[0]
    output_bedpe_file = args[1]
    filter_chimeras(input_bedpe_file, output_bedpe_file)

if __name__ == '__main__': main()