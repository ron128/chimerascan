'''
Created on Jan 29, 2011

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
import os
import subprocess

from find_discordant_reads import DiscordantFragment

SORT_MEMORY = "2g"

def sort_discordant_reads(input_file, output_file):
    outfh = open(output_file, "w")
    ref1_key = "-k%d,%d" % (DiscordantFragment.REF1_COL,
                            DiscordantFragment.REF1_COL)
    ref2_key = "-k%d,%d" % (DiscordantFragment.REF2_COL,
                            DiscordantFragment.REF2_COL)
    tmp_dir = os.path.dirname(output_file)
    args = ["sort", "-S", SORT_MEMORY, "-T", tmp_dir,
            "-t", "\t", "-s", ref1_key, ref2_key,
            input_file]
    logging.debug("sort command args: %s" % args)
    subprocess.call(args, stdout=outfh)
    outfh.close()

def main():
    from optparse import OptionParser
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = OptionParser("usage: %prog <chimeras.bedpe> <chimeras.sorted.bedpe>")
    options, args = parser.parse_args()
    input_file = args[0]
    output_file = args[1]
    logging.info("Sorting discordant BEDPE file")
    logging.info("Input file: %s" % (input_file))
    logging.info("Output file: %s" % (output_file))
    sort_discordant_reads(input_file, output_file)

if __name__ == '__main__':
    main()
