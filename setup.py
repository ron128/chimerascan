'''
chimerascan

Created on Jan 5, 2011

@author: mkiyer
'''
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

import os
import glob

# local imports
import version

# ---- Extension Modules ----------------------------------------------------

def get_pysam_extension_modules():
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
    tabix = Extension("chimerascan.pysam.ctabix", # name of extension
                      ["chimerascan/pysam/ctabix.pyx" ]  +\
                      glob.glob(os.path.join("chimerascan", "pysam", "tabix", "*.c")),
                      library_dirs=[],
                      include_dirs=[ "chimerascan/pysam/tabix", "chimerascan/pysam" ],
                      libraries=[ "z", ],
                      language="c",
                      )
    return samtools, tabix

def get_bx_extension_modules():
    # Interval clustering                
    bx_cluster = Extension("chimerascan.bx.cluster", 
                           ["chimerascan/bx/cluster.pyx", "chimerascan/bx/intervalcluster.c"], 
                           include_dirs=["chimerascan/bx"])
    # Interval intersection
    bx_interval = Extension("chimerascan.bx.intersection",
                            ["chimerascan/bx/intersection.pyx" ])
    return bx_cluster, bx_interval

def get_extension_modules():
    extensions = []
    extensions.extend(get_bx_extension_modules())
    extensions.extend(get_pysam_extension_modules())
    return extensions

def main():
    setup(name="chimerascan",
          version=version.__version__,
          description="chimeric transcript discovery from RNA-seq",
          long_description=__doc__,
          author="Matthew Iyer",
          author_email="mkiyer@umich.edu",
          license="GPL3",
          platforms="ALL",
          url = "http://chimerascan.googlecode.com",
          packages=["chimerascan"],
          ext_modules= get_extension_modules(),
          cmdclass= {'build_ext': build_ext})

if __name__ == '__main__':
    main()