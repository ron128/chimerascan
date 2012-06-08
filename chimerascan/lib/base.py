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
import gzip
import bz2
import zipfile

#
# constants used for library type
#
class LibraryTypes:
    FR_UNSTRANDED = "fr-unstranded"
    FR_FIRSTSTRAND = "fr-firststrand"
    FR_SECONDSTRAND = "fr-secondstrand"

    @staticmethod
    def choices():
        return (LibraryTypes.FR_UNSTRANDED,
                LibraryTypes.FR_FIRSTSTRAND,
                LibraryTypes.FR_SECONDSTRAND)

    @staticmethod
    def same_strand(library_type):
        return (library_type[0] == library_type[1])

imin2 = lambda a,b: a if a <= b else b

def detect_format(f):
    if f.endswith(".gz") or f.endswith(".z"):
        return "gz"
    elif f.endswith(".bz2"):
        return "bz2"
    elif f.endswith(".zip"):
        return "zip"
    else:
        return "txt"

def open_compressed(f):    
    compression_format = detect_format(f)    
    if compression_format == "gz":
        fh = gzip.open(f, "r")
    elif compression_format == "bz2":
        fh = bz2.BZ2File(f, "r")
    elif compression_format == "zip":
        fh = zipfile.ZipFile(f, "r")
    else:
        fh = open(f, "r")
    return fh

def parse_lines(line_iter, numlines=1):
    """
    generator that returns list of 'numlines' lines at a time
    """
    try:
        while True:
            yield [line_iter.next().rstrip() for x in xrange(numlines)]
    except StopIteration:
        pass

def parse_bool(s):    
    return True if s[0].lower() == "t" else False

def parse_string_none(s):
    return None if s == "None" else s

def make_temp(base_dir, suffix=''):
    fd,name = tempfile.mkstemp(suffix=suffix, prefix='tmp', dir=base_dir)
    os.close(fd)
    return name

def check_executable(filename):
    # check that samtools binary exists
    devnullfh = open(os.devnull, 'w')        
    try:
        subprocess.call([filename], stdout=devnullfh, stderr=devnullfh)
    except OSError:
        return False
    devnullfh.close()
    return True

def up_to_date(outfile, infile, nzsize=True):
    if not os.path.exists(infile):
        return False
    if not os.path.exists(outfile):
        return False
    if nzsize and (os.path.getsize(outfile) == 0):
        return False
    return os.path.getmtime(outfile) >= os.path.getmtime(infile)

# in-place XML prettyprint formatter
def indent_xml(elem, level=0):
    i = "\n" + level*"  "
    if len(elem):
        if not elem.text or not elem.text.strip():
            elem.text = i + "  "
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
        for elem in elem:
            indent_xml(elem, level+1)
        if not elem.tail or not elem.tail.strip():
            elem.tail = i
    else:
        if level and (not elem.tail or not elem.tail.strip()):
            elem.tail = i
