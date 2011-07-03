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

def get_header_row():
    return ["5' chrom", "5' breakpoint pos", "5' strand",
            "3' chrom", "3' breakpoint pos", "3' strand",
            "5' transcripts", "3' transcripts",
            "5' genes", "3' genes",
            "Type", "5' -> 3' distance",
            "Multimapping-adjusted Frags",
            "Total Frags", "Spanning Frags",
            "Unique alignment positions",
            "Unique spanning alignment positions"]

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
        for i in xrange(len(fields)):
            fields[i] = ("string", fields[i])
        fields[txs5p_col_num] = ("list", fields[txs5p_col_num][1].split(","))
        fields[txs3p_col_num] = ("list", fields[txs3p_col_num][1].split(","))
        fields[genes5p_col_num] = ("genecards", fields[genes5p_col_num][1].split(","))
        fields[genes3p_col_num] = ("genecards", fields[genes3p_col_num][1].split(","))
        yield fields

def make_html_table(input_file, 
                    show_read_throughs=False):
    line_iter = open(input_file)
    header_line = line_iter.next()[1:]
    header_fields = header_line.strip().split('\t')
    row_iter = generate_row_data(line_iter, 
                                 show_read_throughs=show_read_throughs,
                                 header_fields=header_fields)
    t = env.get_template("table_template.html")
    htmlstring = t.render(colnames=get_header_row(),
                          rows=row_iter)
    return htmlstring

def main():
    from optparse import OptionParser
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = OptionParser("usage: %prog [options] <chimeras.bedpe>")
    parser.add_option("-o", dest="output_file", default=None,
                      help="output file [default=stdout]")
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
    res = make_html_table(input_file, options.show_read_throughs)
    print >>fileh, res
    if options.output_file is not None:
        fileh.close()    


if __name__ == '__main__':
    main()