'''
Created on Jan 5, 2011

@author: mkiyer
'''
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

# ---- Extension Modules ----------------------------------------------------
def get_extension_modules():
    extensions = []
    # Interval clustering                
    extensions.append( Extension( "lib.bx.cluster", [ "lib/bx/cluster.pyx", "lib/bx/intervalcluster.c"], 
                                  include_dirs=["lib/bx"]) )
    # Interval intersection
    extensions.append( Extension( "lib.bx.intersection", [ "lib/bx/intersection.pyx" ] ) )
    return extensions

def main():
    setup(name = "chimerascan",
          ext_modules = get_extension_modules(),
          author = "Matthew Iyer",
          author_email = "mkiyer@umich.edu",
          description = "chimeric transcript discovery tool",
          url = "http://chimerascan.googlecode.com",
          cmdclass= {'build_ext': build_ext})

if __name__ == '__main__': main()