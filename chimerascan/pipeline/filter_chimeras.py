'''
Created on Jan 23, 2011

@author: mkiyer
'''
import logging
import operator
import subprocess
import tempfile
import os

from chimerascan.lib.stats import EmpiricalCdf3D
from nominate_chimeras import Chimera

class SpanningChimera(Chimera):
    def __init__(self):
        Chimera.__init__(self)
        self.spanning_reads = 0
        self.encomp_and_spanning = 0
        self.total_reads = 0
        self.junction_hist = None
        self.spanning_seq = None
        self.spanning_ids = None
        
    def from_list(self, fields):
        # get the chimera fields
        Chimera.from_list(self, fields)
        self.spanning_reads = int(fields[21])
        self.encomp_and_spanning = int(fields[22])
        self.total_reads = int(fields[23])
        self.junction_hist = map(int, fields[24].split(','))
        self.spanning_seq = fields[25]
        self.spanning_ids = fields[26]

    def to_list(self):
        fields = Chimera.to_list(self)
        fields.extend([self.spanning_reads, 
                       self.encomp_and_spanning, 
                       self.total_reads, 
                       ','.join(map(str, self.junction_hist)), 
                       self.spanning_seq, 
                       self.spanning_ids])
        return fields

    @staticmethod
    def parse(line_iter):
        for line in line_iter:
            fields = line.strip().split('\t')
            c = SpanningChimera()
            c.from_list(fields)
            yield c

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
    chimeras = []
    for c in SpanningChimera.parse(open(tmpfile)):
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
    for chimera in chimeras:
        if chimera.distance == 0:
            continue
        yield chimera

def filter_insert_size(chimeras, max_isize):
    for c in chimeras:
        if (c.mate5p.isize + c.mate3p.isize) <= (2*max_isize):
            yield c
        else:
            logging.warning("Removed %s due to insert size %d + %d > %d" %
                            (c.name, c.mate5p.isize, c.mate3p.isize, 2*max_isize))

def get_max_anchor(c, anchor_min):
    # histogram maximum length spanning read
    for maxspan in xrange(len(c.junction_hist)-1, -1, -1):
        if c.junction_hist[maxspan] > 0:
            break
    if maxspan <= anchor_min:
        return 0
    return maxspan

def get_kl_divergence(arr):
    from math import log
    t = sum(arr)
    if t == 0:
        return 0
    expected = t / float(len(arr))
    kldiv = sum((x/float(t))*log(x/expected) for x in arr
                if x > 0)
    return kldiv

def parse_chimera_data(chimeras, anchor_min):
    for c in chimeras:
        maxspan = get_max_anchor(c, anchor_min)
        yield c.encompassing_reads, c.encomp_and_spanning, maxspan
    
def filter_chimeras(input_bedpe_file, 
                    output_bedpe_file,
                    isoforms=False,
                    overlap=False,
                    max_isize=None,
                    anchor_min=4,
                    prob=0.05):    
    # score chimeras
    # build a empirical distribution functions for the chimeras
    logging.info("Determining empirical distribution of chimeras")    
    ecdf = EmpiricalCdf3D(parse_chimera_data(SpanningChimera.parse(open(input_bedpe_file)), 
                                             anchor_min))
    # add filters
    if not isoforms:
        logging.debug("Choosing highest coverage isoforms for each gene pair")        
        iter = filter_isoforms(input_bedpe_file, os.path.dirname(output_bedpe_file))
    else:
        iter = SpanningChimera.parse(open(input_bedpe_file))
    if not overlap:
        logging.debug("Filtering overlapping genes")        
        iter = filter_overlapping_genes(iter)    
    if max_isize is not None:
        logging.debug("Excluding chimeras that violate insert size constraints")        
        iter = filter_insert_size(iter, max_isize=max_isize)
    # get results of filters
    chimeras = list(iter)
    logging.info("Initial filters yielded %d chimeras" % (len(chimeras)))
    # score the chimeras using the distribution function
    logging.info("Scoring chimeras")
    chimera_scores = []
    for c in chimeras:
        maxspan = get_max_anchor(c, anchor_min)
        p = ecdf(c.encompassing_reads, c.encomp_and_spanning, maxspan)
        chimera_scores.append((1.0 - p, c))
    del chimeras
    # sort chimeras
    logging.info("Sorting chimeras by empirical probability")
    chimera_scores = sorted(chimera_scores, key=operator.itemgetter(0))
    # output chimeras
    logging.info("Writing chimeras with p <= %f" % (prob))
    outfh = open(output_bedpe_file, "w")    
    for p,c in chimera_scores:
        #if p <= prob:
        kldiv = get_kl_divergence(c.junction_hist)
        print >>outfh, '\t'.join(['\t'.join(map(str,c.to_list())), str(p), str(kldiv)])
    outfh.close()


def main():
    from optparse import OptionParser
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")    
    parser = OptionParser("usage: %prog [options] <chimeras.bedpe> <filtered_chimeras.bedpe>")
    parser.add_option("--overlap", action="store_true", default=False,
                      help="Retain overlapping genes in the output")
    parser.add_option("--isoforms", action="store_true", default=False,
                      help="Retain all isoforms in the output")    
    parser.add_option("--max-fragment-length", type="int", default=None,
                      help="Maximum allowable insert size")
    parser.add_option("--empirical-prob", type="float", dest="prob", default=0.05,
                      help="Empirical probability of observing reads in data "
                      "(applied after all other filters)")
    parser.add_option("--anchor-min", type="int", default=4,
                      help="Minimum spanning distance to be considered a valid spanning read")
    options, args = parser.parse_args()
    input_bedpe_file = args[0]
    output_bedpe_file = args[1]
    filter_chimeras(input_bedpe_file, 
                    output_bedpe_file, 
                    isoforms=options.isoforms,
                    overlap=options.overlap,
                    max_isize=options.max_fragment_length,
                    anchor_min=options.anchor_min,
                    prob=options.prob)

if __name__ == '__main__': main()

