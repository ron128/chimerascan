'''
chimerascan

Created on Jan 5, 2011

@author: mkiyer
'''
from distutils.core import setup
from distutils.extension import Extension

import os
import glob

# local imports
import version

# ------ Setup instructions -------------------------------------------------

setup_kwargs = {"name": "chimerascan",
                "version": version.__version__,
                "description": "chimeric transcript discovery from RNA-seq",
                "long_description": __doc__,
                "author": "Matthew Iyer",
                "author_email": "mkiyer@umich.edu",
                "license": "GPL3",
                "platforms": "Linux",
                "url": "http://chimerascan.googlecode.com",
                "packages": ["chimerascan",
                             "chimerascan.pysam",
                             "chimerascan.bx",
                             "chimerascan.pipeline",
                             "chimerascan.lib",
                             "chimerascan.tools"],
                "package_data": {'chimerascan.tools': ['table_template.html']},                             
                "scripts": ["chimerascan/chimerascan_run.py",
                            "chimerascan/chimerascan_index.py",
                            "chimerascan/chimerascan_cluster.py",
                            "chimerascan/tools/chimerascan_html_table.py"]}

# ---- Extension Modules ----------------------------------------------------

def get_cython_extension_modules():
    # pysam - samtools
    samtools = Extension("chimerascan.pysam.csamtools", # name of extension
                         ["chimerascan/pysam/csamtools.pyx",
                          "chimerascan/pysam/pysam_util.c"] +\
                          glob.glob( os.path.join( "chimerascan", "pysam", "samtools", "*.c" )),
                          library_dirs=[],
                          include_dirs=[ "chimerascan/pysam/samtools", "chimerascan/pysam" ],
                          libraries=[ "z", ],
                          language="c",
                          define_macros = [('FILE_OFFSET_BITS','64'),
                                           ('_USE_KNETFILE','')])     
    # pysam - tabix
    tabix = Extension("chimerascan.pysam.ctabix", # name of extension
                      ["chimerascan/pysam/ctabix.pyx" ]  +\
                      glob.glob(os.path.join("chimerascan", "pysam", "tabix", "*.c")),
                      library_dirs=[],
                      include_dirs=[ "chimerascan/pysam/tabix", "chimerascan/pysam" ],
                      libraries=[ "z", ],
                      language="c",
                      )
    # Interval clustering                
    bx_cluster = Extension("chimerascan.bx.cluster", 
                           ["chimerascan/bx/cluster.pyx", "chimerascan/bx/intervalcluster.c"], 
                           include_dirs=["chimerascan/bx"])
    # Interval intersection
    bx_interval = Extension("chimerascan.bx.intersection",
                            ["chimerascan/bx/intersection.pyx" ])
    return [samtools, tabix, bx_cluster, bx_interval]

def get_c_extension_modules():
    # pysam - samtools
    samtools = Extension("chimerascan.pysam.csamtools", # name of extension
                         ["chimerascan/pysam/csamtools.c",
                          "chimerascan/pysam/pysam_util.c"] +\
                          glob.glob( os.path.join( "chimerascan", "pysam", "samtools", "*.c" )),
                          library_dirs=[],
                          include_dirs=[ "chimerascan/pysam/samtools", "chimerascan/pysam" ],
                          libraries=[ "z", ],
                          language="c",
                          define_macros = [('FILE_OFFSET_BITS','64'),
                                           ('_USE_KNETFILE','')])     
    # pysam - tabix
    tabix = Extension("chimerascan.pysam.ctabix", # name of extension
                      ["chimerascan/pysam/ctabix.c" ]  +\
                      glob.glob(os.path.join("chimerascan", "pysam", "tabix", "*.c")),
                      library_dirs=[],
                      include_dirs=[ "chimerascan/pysam/tabix", "chimerascan/pysam" ],
                      libraries=[ "z", ],
                      language="c",
                      )
    # Interval clustering                
    bx_cluster = Extension("chimerascan.bx.cluster", 
                           ["chimerascan/bx/cluster.c", "chimerascan/bx/intervalcluster.c"], 
                           include_dirs=["chimerascan/bx"])
    # Interval intersection
    bx_interval = Extension("chimerascan.bx.intersection",
                            ["chimerascan/bx/intersection.c"])
    return [samtools, tabix, bx_cluster, bx_interval]

def main():
    setup(ext_modules=get_c_extension_modules(),
          **setup_kwargs)

if __name__ == '__main__':
    main()