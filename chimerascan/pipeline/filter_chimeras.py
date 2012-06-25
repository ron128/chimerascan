'''
Created on Jan 31, 2011

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
import sys
import os
import logging
import argparse

from chimerascan.lib.chimera import Chimera
from chimerascan.lib import config

def filter_chimeras(input_file, output_file,
                    filter_num_frags,
                    filter_allele_fraction,
                    mask_biotypes,
                    mask_rnames):
    logging.debug("Filtering chimeras")
    logging.debug("\tfragments: %f" % (filter_num_frags))
    logging.debug("\tallele fraction: %f" % (filter_allele_fraction))
    logging.debug("\tmask biotypes: %s" % (','.join(sorted(mask_biotypes))))
    logging.debug("\tmask references: %s" % (','.join(sorted(mask_rnames))))
    # filter chimeras
    num_chimeras = 0
    num_kept_chimeras = 0    
    f = open(output_file, "w")   
    for c in Chimera.parse(open(input_file)):
        num_chimeras += 1
        # number of fragments
        if c.num_frags < filter_num_frags:
            continue
        # allele fraction
        allele_fraction_5p = float(c.num_frags) / (c.num_discordant_frags_5p + c.num_concordant_frags_5p)
        allele_fraction_3p = float(c.num_frags) / (c.num_discordant_frags_3p + c.num_concordant_frags_3p)
        allele_fraction = min(allele_fraction_5p, allele_fraction_3p)
        if allele_fraction < filter_allele_fraction:
            continue
        # masked biotypes and references
        if len(mask_biotypes.intersection(c.biotypes_5p)) > 0:
            continue
        if len(mask_biotypes.intersection(c.biotypes_3p)) > 0:
            continue
        if c.rname5p in mask_rnames:
            continue
        if c.rname3p in mask_rnames:
            continue
        print >>f, str(c)
        num_kept_chimeras += 1
    f.close()
    logging.debug("Total chimeras: %d" % num_chimeras)
    logging.debug("Kept chimeras: %d" % num_kept_chimeras)
    return config.JOB_SUCCESS

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument("--num-frags", dest="num_frags", 
                        type=float, default=config.DEFAULT_FILTER_FRAGS)
    parser.add_argument("--allele-fraction", type=float, 
                        default=config.DEFAULT_FILTER_ALLELE_FRACTION, 
                        dest="allele_fraction", metavar="X",
                        help="Filter chimeras with expression less than "
                        "the specified fraction of the total expression "
                        "level [default=%(default)s")
    parser.add_argument("--mask-biotypes-file", dest="mask_biotypes_file", default=None) 
    parser.add_argument("--mask-rnames-file", dest="mask_rnames_file", default=None)
    parser.add_argument("input_file")
    parser.add_argument("output_file")
    args = parser.parse_args()
    mask_biotypes = set()
    if args.mask_biotypes_file is not None:
        mask_biotypes.update([line.strip() for line in open(args.mask_biotypes_file)])
    mask_rnames = set()
    if args.mask_rnames_file is not None:
        mask_rnames.update([line.strip() for line in open(args.mask_rnames_file)])
    return filter_chimeras(args.input_file, args.output_file,
                           filter_num_frags=args.num_frags,
                           filter_allele_fraction=args.allele_fraction,
                           mask_biotypes=mask_biotypes,
                           mask_rnames=mask_rnames)

if __name__ == "__main__":
    sys.exit(main())