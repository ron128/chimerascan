'''
Created on Feb 12, 2011

@author: mkiyer

chimerascan: chimeric transcript discovery using RNA-seq

Copyright (C) 2011 Matthew Iyer

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''
import logging
import sys
from jinja2 import Environment, PackageLoader

# local imports
from chimerascan.lib.chimera import Chimera, ChimeraTypes

# setup html template environment
env = Environment(loader=PackageLoader("chimerascan", "tools"))

# URLs for special links
GENECARDS_URL = "http://www.genecards.org/cgi-bin/carddisp.pl?gene="
UCSC_POS_URL = "http://genome.ucsc.edu/cgi-bin/hgTracks?"

def get_header_row():
    return ["5' genomic region", "5' strand",
            "3' genomic region", "3' strand",
            "5' transcripts", "3' transcripts",
            "5' genes", "3' genes",
            "Type", "5' -> 3' distance",
            "Multimapping-adjusted Frags",
            "Total Frags", "Spanning Frags",
            "Unique alignment positions",
            "Unique spanning alignment positions",
            "Chimera IDs"]

def generate_row_data(line_iter, show_read_throughs, 
                      header_fields):
    type_col_num = header_fields.index("type") 
    txs5p_col_num = header_fields.index("transcript_ids_5p")
    txs3p_col_num = header_fields.index("transcript_ids_3p")
    genes5p_col_num = header_fields.index("genes5p")
    genes3p_col_num = header_fields.index("genes3p")
    for line in line_iter:
        fields = line.strip().split('\t')
        if ((not show_read_throughs) and 
            (fields[type_col_num] == ChimeraTypes.READTHROUGH)):
            continue
        newfields = []
        # 5' position (chr12:65432) and strand
        newfields.append(("ucsc_pos", ["%s:%s-%s" % (fields[0], fields[1], fields[2])]))
        newfields.append(("string", fields[3]))
        # 3' position (chr12:76543) and strand
        newfields.append(("ucsc_pos", ["%s:%s-%s" % (fields[4], fields[5], fields[6])]))
        newfields.append(("string", fields[7]))
        # transcripts
        newfields.append(("ucsc_pos", fields[txs5p_col_num].split(",")))
        newfields.append(("ucsc_pos", fields[txs3p_col_num].split(",")))
        # genes
        newfields.append(("genecards", fields[genes5p_col_num].split(",")))
        newfields.append(("genecards", fields[genes3p_col_num].split(",")))
        # all but last column not modified
        for i in xrange(genes3p_col_num+1, len(fields)-1):
            newfields.append(("string", fields[i]))
        # last column is a list
        newfields.append(("list", fields[-1].split(",")))
        yield newfields

def make_html_table(input_file, 
                    ucsc_db,
                    show_read_throughs=False):    
    ucsc_pos_url = UCSC_POS_URL + "db=%s&position=" % (ucsc_db)
    line_iter = open(input_file)
    header_line = line_iter.next()[1:]
    header_fields = header_line.strip().split('\t')
    row_iter = generate_row_data(line_iter, 
                                 show_read_throughs=show_read_throughs,
                                 header_fields=header_fields)
    t = env.get_template("table_template.html")
    htmlstring = t.render(colnames=get_header_row(),
                          ucsc_pos_url=ucsc_pos_url,
                          genecards_url=GENECARDS_URL,
                          rows=row_iter)
    return htmlstring

def main():
    from optparse import OptionParser
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = OptionParser("usage: %prog [options] <chimeras.txt>")
    parser.add_option("-o", dest="output_file", default=None,
                      help="output file [default=stdout]")
    parser.add_option("--ucsc-db", dest="ucsc_db", default="hg19",
                      help="UCSC Genome Version (specific to organism and "
                      "revision e.g. 'hg19'")
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
    res = make_html_table(input_file, 
                          ucsc_db=options.ucsc_db,
                          show_read_throughs=options.show_read_throughs)
    print >>fileh, res
    if options.output_file is not None:
        fileh.close()    


if __name__ == '__main__':
    main()