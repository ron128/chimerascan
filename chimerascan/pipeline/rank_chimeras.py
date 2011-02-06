'''
Created on Feb 5, 2011

@author: mkiyer
'''
import logging
import operator
import numpy as np

# local imports
from merge_spanning_alignments import SpanningChimera

PERMISCUITY_THRESHOLD = 0.005

def get_spanning_read_score(c):
    junc_pos = c.junc_pos
    score = 0
    for r in c.spanning_reads:
        anchor = min(junc_pos - r.pos, r.aend - junc_pos)
        score += max(0, anchor) / float(r.mappings)
    return score

def get_ranking_props(c):
    return (c.weighted_cov,
            c.encomp_and_spanning,
            get_spanning_read_score(c),
            int(min(c.mate5p.frac, c.mate3p.frac) > 0.01))       

def hist_interp_prob(H, E, X):
    lo_slices = []
    hi_slices = []
    fracs = []
    for d,edges in enumerate(E):        
        # find bins
        right_ind = np.searchsorted(edges, X[d], side="left")
        left_ind = right_ind - 1
        if right_ind == 0:
            lo_slices.append(slice(0, 0))        
            hi_slices.append(slice(0, 1))
            left_edge = 0        
        else:
            # add all histogram counts up to lo_ind
            lo_slices.append(slice(0, left_ind+1))        
            hi_slices.append(slice(0, right_ind+1))
            left_edge = edges[left_ind]
        # find the fraction between lo_ind and hi_ind
        right_edge = edges[right_ind]
        if right_edge == 0:
            frac = 0
        else:
            frac = (X[d] - left_edge) / (right_edge - left_edge)
        fracs.append(frac)
    # add the initial sum first
    lowval = np.sum(H[lo_slices])
    hival = np.sum(H[hi_slices])
    avgfrac = np.mean(fracs)
    val = lowval + avgfrac * (hival - lowval)
    return val


def rank_chimeras(input_file, output_file, prob):
    '''
    rank the chimeras according to the empirical distribution
    of encompassing read coverage, spanning read coverage, 
    and junction permiscuity
    '''
    # profile the chimeras
    arr = []
    for c in SpanningChimera.parse(open(input_file)):        
        arr.append(get_ranking_props(c))
    arr = np.array(arr)
    # choose bin sizes
    maxbins = 100
    bins = np.zeros(arr.shape[1])
    for d in xrange(arr.shape[1]):
        bins[d] = min(maxbins, len(np.unique(arr[:,d])))
    #H, edges = np.histogramdd(arr, bins=(10,10,10,10), range, normed, weights)
    H, edges = np.histogramdd(arr, bins=bins)
    N = np.sum(H)
    # now rank each chimera using the empirical distribution
    chimera_scores = []
    for c in SpanningChimera.parse(open(input_file)):
        props = get_ranking_props(c)
        p = hist_interp_prob(H, edges, props)
        if p >= prob:
            chimera_scores.append((p/N, c))
    outfh = open(output_file, "w")
    for p,c in sorted(chimera_scores, key=operator.itemgetter(0)):
        print >>outfh, '\t'.join(map(str, c.to_list() + [p]))
    outfh.close() 


def main():
    from optparse import OptionParser
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = OptionParser("usage: %prog [options] <sortedchimeras.bedpe> <chimeras.txt>")
    parser.add_option("--prob", type="float", metavar="p", dest="prob", 
                      default=0.0, help="empirical probability threshold "
                      " for outputting chimeras")
    options, args = parser.parse_args()
    input_file = args[0]
    output_file = args[1]
    rank_chimeras(input_file, output_file, options.prob)

if __name__ == "__main__":
    main()