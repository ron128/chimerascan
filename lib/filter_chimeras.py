'''
Created on Jan 23, 2011

@author: mkiyer
'''
import collections
from collections import defaultdict
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
        c.range5p = map(int, fields[12:14])
        c.range3p = map(int, fields[14:16])
        c.exons5p = fields[16]
        c.exons3p = fields[17]
        c.encompassing_ids = fields[18]
        c.encompassing_seq5p = fields[19]
        c.encompassing_seq3p = fields[20]
        c.spanning_reads = int(fields[21])
        c.encomp_and_spanning = int(fields[22])
        c.total_reads = int(fields[23])
        c.junction_hist = map(int, fields[24].split(','))
        c.spanning_seq = fields[25]
        c.spanning_ids = fields[26]
        return c
    
    def to_tabular(self):
        fields = [self.tx5p, self.start5p, self.end5p,
                  self.tx3p, self.start3p, self.end3p,
                  self.name, self.encompassing_reads,
                  self.strand5p, self.strand3p, self.chimera_type,
                  self.distance, self.range5p[0], self.range5p[1],
                  self.range3p[0], self.range3p[1], 
                  self.exons5p, self.exons3p,
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

def filter_hard_coverage(chimeras, coverage):
    logging.debug("Filtering coverage < %d" % (coverage))
    for chimera in chimeras:
        if chimera.total_reads < coverage:
            continue
        yield chimera

def find_empirical_cutoff(dict, cutoff):
    sorted_keys = sorted(dict)
    total_count = sum(dict[k] for k in sorted_keys)
    print dict
    prob = 0.0
    for bin in sorted_keys:
        if prob >= cutoff:
            break
        count = dict[bin]
        prob += count / float(total_count) 
    return bin, prob
    
def dict_to_probs(self, dict):
    sorted_keys = sorted(dict)
    prob_dict = {}
    n = sum(dict[k] for k in sorted_keys)
    prob = 0.0
    for bin in sorted_keys:
        count = dict[bin]
        prob_dict[bin] = prob
        prob += count / float(n) 
    return n, prob_dict

def profile_empirical_coverage(chimeras, prob=0.95):
    # profile results
    for c in chimeras:
        # histogram encompassing reads per chimera
        
        encomp_dict[c.encompassing_reads] += 1
        encomp_span_dict[c.encomp_and_spanning] += 1
        # histogram maximum length spanning read
        for i in xrange(len(c.junction_hist)-1, -1, -1):
            if c.junction_hist[i] > 0:
                break
        spanning_dict[i] += 1
    # find encompassing and spanning read cutoffs
    encomp_cutoff, encomp_prob = find_empirical_cutoff(encomp_dict, prob)
    span_cutoff, span_prob = find_empirical_cutoff(spanning_dict, prob)
    encomp_span_cutoff, encomp_span_prob = find_empirical_cutoff(encomp_span_dict, prob)
    logging.info("Encompassing reads >= %d has p=%f" % (encomp_cutoff, encomp_prob))
    logging.info("Spanning reads >= %d has p=%f" % (span_cutoff, span_prob))
    logging.info("Encompassing/Spanning reads >= %d has p=%f" % (encomp_span_cutoff, encomp_span_prob))
    return encomp_cutoff, span_cutoff, encomp_span_cutoff

#def profile_empirical_coverage(chimeras, prob=0.95):
#    encomp_dict = collections.defaultdict(lambda: 0)
#    spanning_dict = collections.defaultdict(lambda: 0)
#    encomp_span_dict = collections.defaultdict(lambda: 0)
#    # profile results
#    for c in chimeras:
#        # histogram encompassing reads per chimera
#        encomp_dict[c.encompassing_reads] += 1
#        encomp_span_dict[c.encomp_and_spanning] += 1
#        # histogram maximum length spanning read
#        for i in xrange(len(c.junction_hist)-1, -1, -1):
#            if c.junction_hist[i] > 0:
#                break
#        spanning_dict[i] += 1
#    # find encompassing and spanning read cutoffs
#    encomp_cutoff, encomp_prob = find_empirical_cutoff(encomp_dict, prob)
#    span_cutoff, span_prob = find_empirical_cutoff(spanning_dict, prob)
#    encomp_span_cutoff, encomp_span_prob = find_empirical_cutoff(encomp_span_dict, prob)
#    logging.info("Encompassing reads >= %d has p=%f" % (encomp_cutoff, encomp_prob))
#    logging.info("Spanning reads >= %d has p=%f" % (span_cutoff, span_prob))
#    logging.info("Encompassing/Spanning reads >= %d has p=%f" % (encomp_span_cutoff, encomp_span_prob))
#    return encomp_cutoff, span_cutoff, encomp_span_cutoff

def filter_insert_size(chimeras, 
                       encomp_cutoff,
                       span_cutoff,
                       encomp_span_cutoff,
                       max_isize):
    for c in chimeras:
        if (c.encompassing_reads >= encomp_cutoff or
            c.spanning_reads >= span_cutoff,
            c.encomp_and_spanning >= encomp_span_cutoff):
            yield c
        else:
            range5p = c.range5p[1] - c.range5p[0]
            range3p = c.range3p[1] - c.range3p[0]
            print range5p + range3p, max_isize
            if range5p + range3p <= (2*max_isize):
                yield c

def get_max_anchor(c):
    # histogram maximum length spanning read
    for maxspan in xrange(len(c.junction_hist)-1, -1, -1):
        if c.junction_hist[maxspan] > 0:
            break        
    return maxspan

def parse_chimera_bedpe(filename):
    for line in open(filename):
        yield Chimera.from_tabular(line)

def parse_chimera_data(filename):
    for c in parse_chimera_bedpe(filename):
        maxspan = get_max_anchor(c)
        yield maxspan, c.encompassing_reads, c.encomp_and_spanning
    
def filter_chimeras(input_bedpe_file, 
                    output_bedpe_file,
                    max_isize,
                    prob=0.95):
    # build a empirical distribution functions for the chimeras
    logging.info("Determining empirical distribution of chimeras")
    from empiricalcdf import EmpiricalCdf3D
    ecdf = EmpiricalCdf3D(parse_chimera_data(input_bedpe_file))
    # score the chimeras using the distribution function
    for c in parse_chimera_bedpe(input_bedpe_file):
        maxspan = get_max_anchor(c)
        p = ecdf(maxspan, c.encompassing_reads, c.encomp_and_spanning)
        print '\t'.join(map(str, [p, c.name, c.encompassing_reads, maxspan, c.encomp_and_spanning]))

    return
    
    # profile the distributions of different statistics
    # for assessing chimera quality and choose cutoffs
    # based on empirical probability
    encomp_cutoff, span_cutoff, encomp_span_cutoff = \
        profile_empirical_coverage(parse_chimera_bedpe(input_bedpe_file), prob)
    return
    chimeras = parse_chimera_bedpe(input_bedpe_file)    
    for c in filter_insert_size(chimeras, 
                                encomp_cutoff=encomp_cutoff,
                                span_cutoff=span_cutoff,
                                encomp_span_cutoff=encomp_span_cutoff,
                                max_isize=max_isize):
        pass
        #print c.to_tabular()    
    return

    # setup multiple filters
    it1 = filter_isoforms(input_bedpe_file, os.path.dirname(output_bedpe_file))
    it2 = filter_overlapping_genes(it1)
    chimeras = list(it2)
    sorted_chimeras = sorted(chimeras, key=operator.attrgetter("encomp_and_spanning", "total_reads"), reverse=True)
    outfh = open(output_bedpe_file, "w")
    for c in sorted_chimeras:
        print >>outfh, c.to_tabular()
    outfh.close()
    return
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


def main():
    from optparse import OptionParser
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")    
    parser = OptionParser("usage: %prog [options] <chimeras.bedpe> <filtered_chimeras.bedpe>")
    parser.add_option("--overlap", action="store_true")
    parser.add_option("--prob", type="float", dest="prob", default=0.95)
    options, args = parser.parse_args()
    input_bedpe_file = args[0]
    output_bedpe_file = args[1]
    filter_chimeras(input_bedpe_file, 
                    output_bedpe_file, 
                    max_isize=500,
                    prob=options.prob)

if __name__ == '__main__': main()