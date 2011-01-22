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
    samtools = Extension("lib.pysam.csamtools", # name of extension
                         ["lib/pysam/csamtools.pyx",
                          "lib/pysam/pysam_util.c"] +\
                          glob.glob( os.path.join( "lib", "pysam", "samtools", "*.c" )),
                          library_dirs=[],
                          include_dirs=[ "lib/pysam/samtools", "lib/pysam" ],
                          libraries=[ "z", ],
                          language="c",
                          define_macros = [('FILE_OFFSET_BITS','64'),
                                           ('_USE_KNETFILE','')])     
    tabix = Extension("lib.pysam.ctabix", # name of extension
                      ["lib/pysam/ctabix.pyx" ]  +\
                      glob.glob(os.path.join("lib", "pysam", "tabix", "*.c")),
                      library_dirs=[],
                      include_dirs=[ "lib/pysam/tabix", "lib/pysam" ],
                      libraries=[ "z", ],
                      language="c",
                      )
    return samtools, tabix

def get_bx_extension_modules():
    # Interval clustering                
    bx_cluster = Extension("lib.bx.cluster", 
                           ["lib/bx/cluster.pyx", "lib/bx/intervalcluster.c"], 
                           include_dirs=["lib/bx"])
    # Interval intersection
    bx_interval = Extension("lib.bx.intersection",
                            ["lib/bx/intersection.pyx" ])
    return bx_cluster, bx_interval

def get_extension_modules():
    extensions = []
    extensions.extend(get_bx_extension_modules())
    extensions.extend(get_pysam_extension_modules())
    return extensions

def main():
    setup(name="chimerascan",
          version=version.__version__,
          description="chimeric transcript discovery tool",
          long_description=__doc__,
          author="Matthew Iyer",
          author_email="mkiyer@umich.edu",
          license="MIT",
          platforms="ALL",
          url = "http://chimerascan.googlecode.com",
          py_modules=["lib/pysam/__init__", 
                      "lib/pysam/Pileup", 
                      "lib/pysam/namedtuple",
                      "lib/pysam/version"],
          ext_modules= get_extension_modules(),
          cmdclass= {'build_ext': build_ext})

if __name__ == '__main__':
    main()