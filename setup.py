'''
chimerascan

Created on Jan 5, 2011

@author: mkiyer
'''
from distutils.core import setup
from distutils.extension import Extension

# local imports
import chimerascan

# ------ Setup instructions -------------------------------------------------

setup_kwargs = {"name": "chimerascan",
                "version": chimerascan.__version__,
                "description": "chimeric transcript discovery from RNA-seq",
                "long_description": __doc__,
                "author": "Matthew Iyer",
                "author_email": "mkiyer@umich.edu",
                "license": "GPL3",
                "platforms": "Linux",
                "url": "http://chimerascan.googlecode.com",
                "packages": ["chimerascan",
                             "chimerascan.bx",
                             "chimerascan.pipeline",
                             "chimerascan.lib",
                             "chimerascan.tools"],
                "package_data": {'chimerascan.tools': ['table_template.html']},                             
                "scripts": ["chimerascan/chimerascan_run.py",
                            "chimerascan/chimerascan_index.py",
                            "chimerascan/tools/chimerascan_html_table.py",
                            "chimerascan/tools/make_false_positive_file.py"]}

# ---- Extension Modules ----------------------------------------------------

def get_cython_extension_modules():
    # Interval clustering                
    bx_cluster = Extension("chimerascan.bx.cluster", 
                           ["chimerascan/bx/cluster.pyx", "chimerascan/bx/intervalcluster.c"], 
                           include_dirs=["chimerascan/bx"])
    # Interval intersection
    bx_interval = Extension("chimerascan.bx.intersection",
                            ["chimerascan/bx/intersection.pyx" ])
    return [bx_cluster, bx_interval]

def get_c_extension_modules():
    # Interval clustering                
    bx_cluster = Extension("chimerascan.bx.cluster", 
                           ["chimerascan/bx/cluster.c", "chimerascan/bx/intervalcluster.c"], 
                           include_dirs=["chimerascan/bx"])
    # Interval intersection
    bx_interval = Extension("chimerascan.bx.intersection",
                            ["chimerascan/bx/intersection.c"])
    return [bx_cluster, bx_interval]

def main():
    setup(ext_modules=get_c_extension_modules(),
          **setup_kwargs)

if __name__ == '__main__':
    main()