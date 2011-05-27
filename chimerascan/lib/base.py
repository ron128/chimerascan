'''
Created on Oct 26, 2010

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
import os
import subprocess
import tempfile
import operator

#
# constants used for "strand"
#
POS_STRAND = '+'
NEG_STRAND = '-'
NO_STRAND = '.'

#
# constants used for library type
#
FR_UNSTRANDED = "fr-unstranded"
FR_FIRSTSTRAND = "fr-firststrand"
FR_SECONDSTRAND = "fr-secondstrand"
FF_UNSTRANDED = "ff-unstranded"
FF_FIRSTSTRAND = "ff-firststrand"
FF_SECONDSTRAND = "ff-secondstrand"

LIBRARY_TYPES = (FR_UNSTRANDED,
                 FR_FIRSTSTRAND,
                 FR_SECONDSTRAND,
                 FF_UNSTRANDED,
                 FF_FIRSTSTRAND,
                 FF_SECONDSTRAND) 

def cmp_strand(a, b):
    '''return True if strands are compatible, False otherwise'''
    if (a == NO_STRAND) or (b == NO_STRAND):
        return True
    return a == b

def up_to_date(outfile, infile):
    if not os.path.exists(infile):
        return False
    if not os.path.exists(outfile):
        return False
    if os.path.getsize(outfile) == 0:
        return False    
    return os.path.getmtime(outfile) >= os.path.getmtime(infile)

def parse_bool(s):    
    return True if s[0].lower() == "t" else False

def parse_string_none(s):
    return None if s == "None" else s

def make_temp(base_dir, suffix=''):
    fd,name = tempfile.mkstemp(suffix=suffix, prefix='tmp', dir=base_dir)
    os.close(fd)
    return name

def get_read_length(fastq_file):
    f = open(fastq_file)
    f.next()
    seq = f.next().strip()
    f.close()
    return len(seq)

def get_read_length_compressed(input_file):
    import gzip
    import bz2    
    suffix = os.path.splitext(input_file)[-1]
    if suffix == '.gz':
        f = gzip.GzipFile(input_file, 'r')
    elif suffix == '.bz2':
        f = bz2.BZ2File(input_file, 'r')
    else:
        f = open(input_file, 'r')
    f.next()
    seq = f.next().strip()
    f.close()
    return len(seq)

def check_executable(filename):
    # check that samtools binary exists
    devnullfh = open(os.devnull, 'w')        
    try:
        subprocess.call([filename], stdout=devnullfh, stderr=devnullfh)
    except OSError:
        return False
    devnullfh.close()
    return True

def get_refs_from_bowtie_index(bowtie_index, split=True):
    args = ['bowtie-inspect', '-s', bowtie_index]    
    p = subprocess.Popen(args, stdout=subprocess.PIPE)
    output = p.communicate()[0]
    refs = []
    for line in output.split('\n'):
        if not line:
            continue
        fields = line.split('\t')
        if fields[0].startswith('Sequence'):
            refname = fields[1].split()[0] if split else fields[1]
            refs.append((refname, int(fields[2])))
    return refs




