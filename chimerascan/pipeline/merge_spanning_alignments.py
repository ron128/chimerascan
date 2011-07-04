'''
Created on Nov 7, 2010

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
import collections

# local imports
from chimerascan import pysam
from chimerascan.lib.chimera import Chimera, DiscordantRead, \
    DiscordantTags, DISCORDANT_TAG_NAME, \
    OrientationTags, ORIENTATION_TAG_NAME
from chimerascan.lib.breakpoint import Breakpoint
from chimerascan.lib.sam import parse_reads_by_qname, parse_pe_reads
from chimerascan.lib.config import GENE_REF_PREFIX

def get_mismatch_positions(md):
    x = 0
    pos = []
    for y in xrange(len(md)):
        if md[y].isalpha():
            offset = int(md[x:y])
            pos.append(offset)
            x = y + 1
    return pos

def check_breakpoint_alignment(r, b,
                               homology5p,
                               homology3p,
                               anchor_min,
                               anchor_length,
                               anchor_mismatches):
    """
    returns True if read 'r' meets criteria for a valid
    breakpoint spanning read, False otherwise
    
    r - pysam AlignedRead object
    b - Breakpoint object
    homology5p - number of bases of 5' homology
    homology3p - number of bases of 3' homology
    """
    # check if read spans breakpoint
    if not (r.pos < b.pos < r.aend):
        return False   
    # calculate amount in bp that read overlaps breakpoint
    # and ensure overlap is sufficient
    left_anchor_bp = b.pos - r.pos
    if left_anchor_bp < max(homology5p, anchor_min):
        return False
    right_anchor_bp = r.aend - b.pos
    if right_anchor_bp < max(homology3p, anchor_min):
        return False
    # ensure that alignments with anchor overlap less than 'anchor_length'
    # do not have more than 'anchor_mismatches' mismatches in the 
    # first 'anchor_length' bases
    if min(left_anchor_bp, right_anchor_bp) < anchor_length:
        # find interval of smallest anchor
        if left_anchor_bp < anchor_length:
            anchor_interval = (0, left_anchor_bp)
        else:
            aligned_length = r.aend - r.pos
            anchor_interval = (aligned_length - right_anchor_bp, aligned_length)      
        # get positions where mismatches occur
        mmpos = get_mismatch_positions(r.opt('MD'))
        # see if any mismatches lie in anchor interval
        anchor_mm = [pos for pos in mmpos
                     if anchor_interval[0] <= pos < anchor_interval[1]]
        if len(anchor_mm) > anchor_mismatches:
            # too many mismatches within anchor position
            return False
    return True

def filter_spanning_reads(reads,
                          tid_breakpoint_dict,
                          homology_dict,
                          anchor_min,
                          anchor_length,
                          anchor_mismatches):
    for i,r in enumerate(reads):
        if r.is_unmapped:
            continue
        # get breakpoint information
        b = tid_breakpoint_dict[r.rname]
        # add tags to read
        r.tags = r.tags + [("HI",i),
                           ("IH",len(reads)),
                           ("NH", len(reads)),
                           (DISCORDANT_TAG_NAME, DiscordantTags.DISCORDANT_GENE),
                           (ORIENTATION_TAG_NAME, OrientationTags.NONE)]
        # determine whether this is breakpoint
        # alignment meets filtering criteria
        for chimera_name in b.chimera_names:
            homology5p, homology3p = homology_dict[chimera_name]
            if check_breakpoint_alignment(r, b,
                                          homology5p,
                                          homology3p,
                                          anchor_min,
                                          anchor_length,
                                          anchor_mismatches):
                yield r,b,chimera_name

def make_tid_breakpoint_dict(bamfh, breakpoint_map_file):
    rname_tid_map = dict((rname,i) for i,rname in enumerate(bamfh.references))
    # make a dictionary to lookup breakpoint information from reference 'tid'
    tid_breakpoint_dict = {}
    for line in open(breakpoint_map_file):
        fields = line.strip().split('\t')
        b = Breakpoint.from_list(fields)
        # add to dictionary
        tid = rname_tid_map[b.name]
        tid_breakpoint_dict[tid] = b
    return tid_breakpoint_dict

def process_encomp_spanning_reads(bam_file,
                                  breakpoint_map_file,
                                  homology_dict,                           
                                  anchor_min, 
                                  anchor_length,
                                  anchor_mismatches):
    """
    returns a dictionary keyed by a unique chimera name (string) with
    and value is a set of read names corresponding to spanning reads
    """
    # map reference names to numeric ids
    bamfh = pysam.Samfile(bam_file, "rb")
    tid_breakpoint_dict = make_tid_breakpoint_dict(bamfh, breakpoint_map_file)
    # assign read alignments to breakpoints
    num_multimaps = 0    
    num_reads = 0
    num_alignments = 0
    num_filtered_hits = 0
    alignment_dict = collections.defaultdict(lambda: [])    
    for reads in parse_reads_by_qname(bamfh):
        # track basic statistics
        num_reads += 1
        if len(reads) > 1:
            num_multimaps += 1
        num_alignments += len(reads)
        # iterate through reads
        for r,b,chimera_name in filter_spanning_reads(reads,
                                                      tid_breakpoint_dict,
                                                      homology_dict,
                                                      anchor_min,
                                                      anchor_length,
                                                      anchor_mismatches):
            dr = DiscordantRead.from_read(r)
            dr.is_spanning = True
            alignment_dict[chimera_name].append(dr)
            num_filtered_hits += 1
    bamfh.close()
    # report statistics
    logging.debug("Encompassing/Spanning Fragments: %d" % (num_reads))
    logging.debug("\tMultimapping: %d" % (num_multimaps))
    logging.debug("\tAlignments: %d" % (num_alignments))
    logging.debug("\tBreakpoint spanning alignments: %d" % (num_filtered_hits))    
    # return dictionary keyed by chimera name with all valid breakpoint
    # reads mapping to that chimera
    return alignment_dict

def parse_unaligned_spanning_bams(unaligned_bamfh,
                                  unaligned_spanning_bamfh):
    unaligned_bam_iter = iter(parse_pe_reads(unaligned_bamfh))
    try:
        pe_reads = unaligned_bam_iter.next()        
    except StopIteration:
        return
    for spanning_reads in parse_reads_by_qname(unaligned_spanning_bamfh):
        # sync the spanning BAM file and the original unaligned BAM file
        while spanning_reads[0].qname != pe_reads[0][0].qname:
            pe_reads = unaligned_bam_iter.next()
        yield pe_reads,spanning_reads

def process_unaligned_spanning_reads(unaligned_bam_file,
                                     spanning_bam_file,
                                     breakpoint_map_file,
                                     homology_dict,
                                     rname_chimeras_5p,
                                     rname_chimeras_3p,                       
                                     anchor_min, 
                                     anchor_length,
                                     anchor_mismatches):   
    unaligned_bamfh = pysam.Samfile(unaligned_bam_file, "rb")
    spanning_bamfh = pysam.Samfile(spanning_bam_file, "rb")
    tid_breakpoint_dict = make_tid_breakpoint_dict(spanning_bamfh, breakpoint_map_file)    
    alignment_dict = collections.defaultdict(lambda: [])    
    num_frags = 0
    num_alignments = 0
    num_filtered_hits = 0    
    for pe_reads,spanning_reads in parse_unaligned_spanning_bams(unaligned_bamfh, spanning_bamfh):
        num_frags += 1
        num_alignments += len(spanning_reads)
        # find which of the original reads was unmapped        
        r1_unmapped = any(r.is_unmapped for r in pe_reads[0])
        r2_unmapped = any(r.is_unmapped for r in pe_reads[1])
        # setup filtering iterator through spanning reads
        filter_iter = filter_spanning_reads(spanning_reads,
                                            tid_breakpoint_dict,
                                            homology_dict,
                                            anchor_min,
                                            anchor_length,
                                            anchor_mismatches)
        if r1_unmapped and r2_unmapped:
            # first possibility is that both reads are unmapped.  this can
            # happen when there is a small fragment with overlapping reads            
            # for this to be a valid spanning fragment, both of the reads 
            # must map to the same breakpoint
            breakpoint_read_dict = collections.defaultdict(lambda: [[],[]])
            # iterate through reads
            for r,b,chimera_name in filter_iter:
                readnum = 0 if r.is_read1 else 1
                breakpoint_read_dict[b.name][readnum].append((r,chimera_name))
            # check for breakpoints that have both mates mapped
            for b_name,b_pe_reads in breakpoint_read_dict.iteritems():
                if (len(b_pe_reads[0]) == 0) or (len(b_pe_reads[1]) == 0):
                    continue
                # this is a good read pair                    
                for readnum,read_chimera_tuples in enumerate(b_pe_reads):
                    for r,chimera_name in read_chimera_tuples:                            
                        # make discordant read object and add to list
                        dr = DiscordantRead.from_read(r)
                        dr.is_spanning = True
                        alignment_dict[chimera_name].append(dr)
                        num_filtered_hits += 1
        else:
            # one of the two reads is unmapped and the other is mapped.  
            # check that the mapped read in the pair aligns to a predicted 
            # chimera
            # get chimeras corresponding to mapped read
            mapped_readind = 0 if r2_unmapped else 1
            mapped_chimeras = set()
            for r in pe_reads[mapped_readind]:
                rname = unaligned_bamfh.references[r.rname]
                is_sense = not r.is_reverse
                if is_sense and (rname in rname_chimeras_5p):
                    mapped_chimeras.update(rname_chimeras_5p[rname])
                elif (not is_sense) and (rname in rname_chimeras_3p):
                    mapped_chimeras.update(rname_chimeras_3p[rname])
            # intersect with chimeras corresponding to breakpoint
            # spanning reads
            for r,b,chimera_name in filter_iter:
                if chimera_name in mapped_chimeras:
                    # make discordant read object and add to list
                    dr = DiscordantRead.from_read(r)
                    dr.is_spanning = True
                    alignment_dict[chimera_name].append(dr)
                    num_filtered_hits += 1
    unaligned_bamfh.close()
    spanning_bamfh.close()
    # report statistics
    logging.debug("Unaligned Spanning Fragments: %d" % (num_frags))
    logging.debug("\tAlignments: %d" % (num_alignments))
    logging.debug("\tBreakpoint spanning alignments: %d" % (num_filtered_hits))    
    return alignment_dict


def merge_spanning_alignments(input_chimera_file, 
                              encomp_spanning_bam_file,
                              unaligned_bam_file,
                              spanning_bam_file, 
                              breakpoint_map_file, 
                              output_chimera_file,
                              anchor_min, 
                              anchor_length,
                              anchor_mismatches):
    # read breakpoint homology information used to
    # filter spanning reads and build a lookup table
    # from reference name to chimera name
    homology_dict = {}
    rname_chimeras_5p = collections.defaultdict(lambda: set())
    rname_chimeras_3p = collections.defaultdict(lambda: set())
    for c in Chimera.parse(open(input_chimera_file)):
        homology_dict[c.name] = (c.breakpoint_homology_5p, 
                                 c.breakpoint_homology_3p)
        rname_chimeras_5p[GENE_REF_PREFIX + c.partner5p.tx_name].add(c.name)
        rname_chimeras_3p[GENE_REF_PREFIX + c.partner3p.tx_name].add(c.name)
    # aggregate breakpoint alignment data from bam file
    # into a dictionary keyed by chimera name
    encomp_spanning_dict = \
        process_encomp_spanning_reads(encomp_spanning_bam_file,
                                      breakpoint_map_file,
                                      homology_dict,
                                      anchor_min, 
                                      anchor_length,
                                      anchor_mismatches)
    unaligned_spanning_dict = \
        process_unaligned_spanning_reads(unaligned_bam_file,
                                         spanning_bam_file,
                                         breakpoint_map_file,
                                         homology_dict,
                                         rname_chimeras_5p,
                                         rname_chimeras_3p,                       
                                         anchor_min, 
                                         anchor_length,
                                         anchor_mismatches)            
    # add breakpoint spanning reads to chimera objects
    logging.debug("Writing spanning hits into chimeras")
    f = open(output_chimera_file, "w")
    for c in Chimera.parse(open(input_chimera_file)):
        # index encomp reads by qname
        encomp_qname_pair_dict = dict((pair[0].qname,pair) for pair in c.encomp_read_pairs)
        spanning_frags = collections.defaultdict(lambda: ([],[]))
        # does chimera have any encompassing/spanning reads?
        if c.name in encomp_spanning_dict:
            for dr in encomp_spanning_dict[c.name]:
                if dr.qname not in encomp_qname_pair_dict:
                    continue
                encomp_qname_pair_dict[dr.qname][dr.readnum].is_spanning = True
                spanning_frags[dr.qname][dr.readnum].append(dr)
                c.spanning_reads.append(dr)                
        # does chimera have unaligned spanning reads?
        if c.name in unaligned_spanning_dict:
            for dr in unaligned_spanning_dict[c.name]:
                spanning_frags[dr.qname][dr.readnum].append(dr)
                c.spanning_reads.append(dr)     
        c.num_spanning_frags = len(spanning_frags.keys())
        # write to output file
        fields = c.to_list()
        print >>f, '\t'.join(map(str, fields))                
    f.close()

def main():
    from optparse import OptionParser
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")    
    parser = OptionParser("usage: %prog [options] <chimeras.in.txt> "
                          "<encomp.bam> <unaligned.bam> <breakpoint_map> "
                          "<chimeras.out.txt>")
    parser.add_option("--anchor-min", type="int", dest="anchor_min", default=4)
    parser.add_option("--anchor-length", type="int", dest="anchor_length", default=8)
    parser.add_option("--anchor-mismatches", type="int", dest="anchor_mismatches", default=0)
    options, args = parser.parse_args()
    input_chimera_file = args[0]
    encomp_spanning_bam_file = args[1]
    unaligned_bam_file = args[2]
    spanning_bam_file = args[3]
    breakpoint_map_file = args[4]
    output_chimera_file = args[5]
    merge_spanning_alignments(input_chimera_file, 
                              encomp_spanning_bam_file,
                              unaligned_bam_file,
                              spanning_bam_file, 
                              breakpoint_map_file, 
                              output_chimera_file,
                              options.anchor_min, 
                              options.anchor_length,
                              options.anchor_mismatches)

if __name__ == '__main__':
    main()
