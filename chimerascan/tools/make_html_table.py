'''
Created on Feb 12, 2011

@author: mkiyer
'''
import logging
import os
import sys
import jinja2
from jinja2 import Environment, FileSystemLoader

# local imports
from chimerascan.pipeline.merge_spanning_alignments import SpanningChimera
from chimerascan.pipeline.nominate_chimeras import CHIMERA_READTHROUGH

env = Environment(loader=FileSystemLoader(os.path.dirname(__file__)))


def get_header_row():
    return ["5' ucsc id",
            "5' exons",
            "3' ucsc id",
            "3' exons",
            "Name",
            "Type",
            "5'-> 3' distance",
            "Empirical pvalue",
            "Weighted coverage",
            "Encompassing reads",
            "Spanning reads",
            "Encompassing and spanning",
            "Total reads",
            "Breakpoint spanning histogram",
            "Multimapping histogram",
            "Breakpoint ID",
            "5'->3' permiscuity",
            "3'->5' permiscuity"]

def generate_row_data(chimeras, prob_cutoff, show_read_throughs):                      
    for c in chimeras:
        prob = float(c.extra_fields[-1])
        if prob > prob_cutoff:
            break
        if ((not show_read_throughs) and 
            (c.chimera_type == CHIMERA_READTHROUGH)):
            continue
        yield [c.mate5p.tx_name,
               '%d-%d' % (c.mate5p.exon_start_num, c.mate5p.exon_end_num),
               c.mate3p.tx_name,
               '%d-%d' % (c.mate3p.exon_start_num, c.mate3p.exon_end_num),
               c.name,
               c.chimera_type,
               c.distance,
               c.extra_fields[-1], # empirical pvalue
               c.weighted_cov,
               c.encompassing_reads,
               c.num_spanning_reads,
               c.encomp_and_spanning,
               c.encomp_or_spanning,
               c.extra_fields[-2], # breakpoint hist
               ','.join(map(str,c.multimap_cov_hist)),
               c.junc_name,
               c.mate5p.frac,
               c.mate3p.frac]

def make_html_table(input_file, 
                    prob_cutoff=0.6,
                    show_read_throughs=False):
    row_iter = generate_row_data(SpanningChimera.parse(open(input_file)), 
                                 prob_cutoff=prob_cutoff,
                                 show_read_throughs=show_read_throughs)
    t = env.get_template("table_template.html")
    htmlstring = t.render(colnames=get_header_row(),
                          rows=row_iter)
    return htmlstring

def main():
    from optparse import OptionParser
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = OptionParser("usage: %prog [options] <chimeras.ranked.bedpe>")
    parser.add_option("-o", dest="output_file", default=None,
                      help="output file [default=stdout]")
    parser.add_option("--empirical-prob", dest="empirical_prob", 
                      type="float", default=0.6,
                      help="probability threshold (0-1) [default=%default]")
    parser.add_option("--read-throughs", dest="show_read_throughs",
                      action="store_true", default=False,
                      help="include read-through chimeras in output "
                      "[default=%default]")
    options, args = parser.parse_args()
    input_file = args[0]
    if options.output_file is None:
        fileh = sys.stdout
    else:
        fileh = open(options.output_file, "w")
    res = make_html_table(input_file, options.empirical_prob,
                          options.show_read_throughs)
    print >>fileh, res
    #rank_chimeras(input_file, output_file, options.empirical_prob)
    if options.output_file is not None:
        fileh.close()    


if __name__ == '__main__':
    main()